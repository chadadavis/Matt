/* pdb.h -- Contains prototypes for loading and manipulating pdb files.
 * Can keep track of most information in a PDB file, recalculate secondary
 * structure (Alpha helices and beta strands only), and has some provisions
 * for handling errors/inconsistencies somewhat intelligently.
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

#ifndef PDB_H
#define PDB_H

#include "FileReader.h"
#include "chain.h"



struct PDB_MODEL {
	int numChains;
	PDBChain **chains;
	char *pdbName;

	int number;
};
typedef struct PDB_MODEL PDBModel;

/* name is pdb file name. */
PDBModel *CreateModel(char *name, int number);
/* Removes a chain from the structure and returns that chain.
 * This way, can pass chains around independent of the rest of
 * the structure, and delete everything else.
 */
PDBChain *PopChain(PDBModel *model, char chain);

PDBChain *GetChain(PDBModel *model, char chain);
/* Just like GetChain, but creates a new chain if it doesn't already exist. */
PDBChain *AddNewChain(PDBModel *model, char chain);

void AddExistingChain(PDBModel *model, PDBChain *chain);

void RotateModel(PDBModel *model, const Matrix *m);
void CleanupModel(PDBModel *model);



struct PDB_BORING_LINES {
	int type;
	char *value;
};
typedef struct PDB_BORING_LINES PDBBoringLines;

struct PDB_DATA {
	PDBBoringLines boringLines;
	char *header;
	char *cryst;
	PDBModel** models;
	int numModels;
	char *name;
};
typedef struct PDB_DATA PDBData;

void RotatePDB(PDBData *pdb, const Matrix *m);
void CleanupPDB(PDBData *pdb);

PDBData *LoadPDBFile(char *name, int seqres, int maxGaps, int caonly, int het);

void AddChain(PDBData *pdb, int model, PDBChain *chain);

void Output(PDBData *pdb, FILE *out, int seqres, int het, int secondary, int extraInfo);
void OutputChains(PDBChain **chain, int numChains, FILE *out, int seqres, int het, int secondary, int extraInfo);
PDBModel *FindModel(PDBData *pdb, int number);
void OutputFasta(PDBData *pdb, FILE *out, int het);
void CalcPDBSecondary(PDBData *pdb);
PDBModel *AddModel(PDBData *pdb, int number);

#endif
