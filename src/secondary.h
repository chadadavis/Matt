/* Protein.h -- Structures that contain secondary structure info, as well as two
 * database structures, to get them out of the way.
 * Two different representations for beta-sheets.  One containing only
 * beta-strand pairs, and one corresponding to PDB beta-sheet entries.
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

#ifndef SECONDARY_H
#define SECONDARY_H
enum HelixType {
	RIGHT_ALPHA = 1,
	RIGHT_OMEGA = 2,
	RIGHT_PI    = 3,
	RIGHT_GAMMA = 4,
	RIGHT_310   = 5,
	LEFT_ALPHA  = 6,
	LEFT_OMEGA  = 7,
	LEFT_GAMMA  = 8,
	HELIX_27    = 9,
	POLYPROLENE = 10
};

#define TYPE_HELIX 2
#define TYPE_SHEET 1

struct ALPHA_HELIX {
	int start;
	int length;
	int type;
	char id[4];
	char comment[31];
};
typedef struct ALPHA_HELIX AlphaHelix;


struct BETA_RESIDUE_PAIR {
	int res1;
	int res2;
	int dir;
};
typedef struct BETA_RESIDUE_PAIR BetaResiduePair;

struct HYDROGEN_BOND {
	int res1, res2;
	int atom1, atom2;
};

typedef struct HYDROGEN_BOND HydrogenBond;

/* Represents an interacting pair of beta strands.
 * No effort is made to assemble them into sheets here,
 * as sheets can split, etc.
 */
struct BETA_PAIR {
	int start1, start2;
	int direction;
	int length;
};
typedef struct BETA_PAIR BetaPair;

int FindFirstMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair);
int FindSecondMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair);
int FindMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair);
void FlipPair(BetaPair *betaPair);

/* Used in beta sheets corresponding to sheets in PDB file */
struct BETA_SHEET_PAIR {
	int res1, res2;
	int length;
	int direction;
	/* Used temporarily, when loading. */
	char insert1, insert2;
	/* Used for spanning chains.
	 * Currently only used when loaded from PDB files.
	 */
	char seq1, seq2;
};
typedef struct BETA_SHEET_PAIR BetaSheetPair;

/* Correspond to sheets in PDB file */
struct BETA_SHEET {
	char id[4];
	BetaSheetPair *strands;
	int numStrands;
};
typedef struct BETA_SHEET BetaSheet;

void FlipSheet(BetaSheet *b);


struct SSBOND {
	int res1, res2;
	int sym1, sym2;
};
typedef struct SSBOND SSbond;

struct CIS_PEP {
	int res1;
	int res2;
	double angle;
};
typedef struct CIS_PEP CisPep;

struct DB_REF {
	int start;
	int end;
	char idCode[5];
	/* everything else.  No need to parse it. */
	char stuff[44];
};
typedef struct DB_REF DBRef;

struct SEQ_ADV {
	int res;
	char idCode[5];
	char resName[4];
	/* everything else.  No need to parse it. */
	char stuff[47];
};
typedef struct SEQ_ADV SeqAdv;
#endif
