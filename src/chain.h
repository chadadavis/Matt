/* chain.h -- Contains prototypes for creating and manipulating pdb chains.
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

#ifndef CHAIN_H
#define CHAIN_H

#include "Vector.h"
#include "secondary.h"

#include "Protein.h"

#include <stdio.h>

struct RESIDUE {
	/* 3-character string containing "residue" name.  Spaces matter, ends with a null. */
	char resName[4];

	/* indices of first and last atoms in this residue. */
	int firstAtom;
	int lastAtom;

	/* resNum and resLetter are the PDB number and letter used to identify the
	   residue.  No letter is kept as a space */
	int resNum;
	char resLetter;
};

typedef struct RESIDUE Residue;

void InitResidue(Residue *res, char *name);


struct ATOM {
	Vector pos;

	float occupancy;
	float tempFactor;

	/* Index of atom - should be unique within chain. */
	int atomIndex;

	/* Number of the residue or het block this is a part of. */
	int res;

	/* Real 0-based residue index. */
	int resIndex;

	/* name of atom. */
	char name[5];

	/* Leading space removed, if any.  All caps. */
	char element[3];

	/* Name of residue from PDB file.  Theoretically
	 * redundant, but sometimes residues and funky het
	 * blocks share a number.
	 */
	char resName[4];

	/* Better known as insertion code. */
	char resLetter;

	char charge;

	char altLoc;
};
typedef struct ATOM Atom;


struct TEMP_RESIDUE {
	int finalIndex;

	int res;
	char resLetter;

	/* if this is in a het block or atom entry. */
	char het;
	/* If this has an alpha carbon. */
	char alphaCarbon;

	/* Name of residue from PDB file.  Temp variable
	   used until atoms aligned with seqres entries. */
	char resName[4];

};
typedef struct TEMP_RESIDUE TempResidue;

/* Note:  For simplicity, chains without a chain id are left with
 * a chain id of space.  Elsewhere, 0 is used instead.
 * The the model commands to find chains will return chains with names
 * of space when a chain id of binary 0 is given.
 */
struct PDB_CHAIN {
	CisPep *cisPeps;
	int numCisPeps;
	SeqAdv *seqAdvs;
	int numSeqAdvs;
	DBRef *dbrefs;
	int numDbrefs;
	SSbond *ssbonds;
	int numSSbonds;

	BetaSheet *betaSheets;
	AlphaHelix *alphaHelices;
	int numBetaSheets;
	int numAlphaHelices;

	HydrogenBond *hbonds;
	BetaPair *betaPairs;
	int numHbonds;
	int numBetaPairs;

	/* Single character storage of sequence. */
	Sequence *seq;
	/* Full listing of all the residues along with numbers of all atoms in each residue. */
	Residue *residues;

	int length;

	Atom *atoms;
	int numAtoms;
	/* Letter identifying this chain.  0 is used if there is no chain.
	  Used for outputting the file. */
	char chainName;

	/* Used to identify the chain in text files.  Contains the file name, less .pdb or .ent, less a leading pdb,
	 * and the chain name, if relevant.
	 */
	char *idString;

	/* variables just used for loading. */
	int terminated;
	TempResidue *tempResidues;
	int numTempResidues;


	/* 1 if we've found any atom records.  Used for guessing if we have a chain
	 * or just a bunch of molecules.  Have a chain when a chain has a TER
	 * record, a SEQRES record or at least one ATOM record.  The last is because
	 * a lot of things that output coordinates only give atom records.  The first two
	 * aren't really needed, but reduce the chance of problems.
	 */
	char tempAtoms;

	int secondaryCalculated;
};
typedef struct PDB_CHAIN PDBChain;

void AddHbond(PDBChain *chain, int res1, int res2);
void AddStrandPair(PDBChain *chain, int start1, int start2, int direction, int length);

int GetAtomIndex(PDBChain *chain, int pos, char*desc);

void CalcBetaSheets(PDBChain *chain);
void CalcStandardHelixType(PDBChain *chain, int skip, int name);
void CalcHelices(PDBChain *chain);
void CalcHbonds(PDBChain *chain);

void TerminateChainOutput(PDBChain *chain, FILE *out, int *index, int res, char c);

/* Note:  None of the following queries will cause problems if a residue
 * index goes beyond the end of the chain in either direction.
 */

/* Returns 1 if the NH in res1 is hydrogen bonded to the O in res2.
 * Backbone residues only.  Otherwise, returns 0
 */
int AreHbondedAsymmetric(PDBChain *chain, int res1, int res2, HydrogenBond **hbond);
/* Returns 1 if the NH or O in res1 is hydrogen bonded to the O or NH in res2.
 * Backbone residues only.  Otherwise, returns 0
 */
int AreHbonded(PDBChain *chain, int res1, int res2, HydrogenBond **hbond);

/* Checks if res1 and res2 are aligned in a beta sheet.  Returns 0 if they aren't,
 * 1 if they are and are parallel, -1 for antiparallel.
 */
int AreAlignedPairs(PDBChain *chain, int res1, int res2);
void CalcChainSecondary(PDBChain *chain);

PDBChain *CreateChain(char *name, char chain);

void RotateChain(PDBChain *chain, const Matrix *m);

void Renumber(PDBChain *chain);

void CleanupChain(PDBChain *chain);

Residue * GetResidue(PDBChain *chain, int pos);
Atom *GetAtom(PDBChain *chain, int pos, char*desc);

/* Adds a residue to the chain, both to residues and to seq. */
void AddResidue(PDBChain *chain, char *resName);

/* Ter is 1 if we run into a ter record, 0 otherwise. */
void Terminate(PDBChain *chain, int maxGaps, int ter);
void End(PDBChain *chain, int maxGaps);

int OutputChainSeqres(PDBChain *chain, FILE *out, int het);
void OutputChainSheets(PDBChain *chain, FILE *out, int het, int *numSheets, char **sheetIDs, int *sheetLines);
void OutputChainAlphas(PDBChain *chain, FILE *out, int het, int *numHelices, char **helixIDs);
void OutputChainSSbonds(PDBChain *chain, FILE *out, int het, int *numSSbonds);
void OutputChainCisPeps(PDBChain *chain, FILE *out, int het, int *numCisPeps, int model);
int OutputChain(PDBChain *chain, FILE *out, int seqres, int het, int *index);
int OutputChainFasta(PDBChain *chain, FILE *out, int het);

/* Converts from pdb position code to position in
 * sequence, used in the sequence and residue arrays.
 * returns 0 on failure, 1 on success.  Result is output
 * to first argument.  Insertion code is optional, and assumed
 * to be the first one, if any codes are in use.  Note that residues
 * given a number by the algorithm have insertion codes of 0, while
 * atoms without a code but with residues in the original file have
 * codes of ' '.  If numbering is not sequential, not guaranteed
 * to give the first occurence.  Takes O(log n) if numbering is sequential
 * and finds something.  Otherwise, takes O(n).
 */
int GetResidueIndex(PDBChain *chain, int *position, int pdbPosition, char insertionCode);

/* Guaranteed to give the first occurance with the given codes. */
int GetResidueFirstIndex(PDBChain *chain, int *position, int pdbPosition, char insertionCode);
/* Does the same thing, but works with a portion of a string from a raw PDB line. */
int GetResidueFirstIndexString(PDBChain *chain, int *position, char * rawString);
void AddAlphaHelix(PDBChain *chain, AlphaHelix *helix);

void RelabelChains(PDBChain **chains, int numChains);
void AssembleSheets(PDBChain *chain);

PDBChain *DuplicateChain(PDBChain *chain);

#endif
