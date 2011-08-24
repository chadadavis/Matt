/* Extend.h -- Contains prototypes to extend an alignment.
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

#ifndef EXTEND_H
#define EXTEND_H
#include "MultipleAlignment.h"

/*  Standard extension phase.  Extends existing blocks, allows up to
 *  5 residue overlap.
 */
void Extend(MultipleAlignment *ma, Matrix **matrices, int useSkips);

/* Final extension pass. Extends using specified rmsd. Doesn't handle
 * missing residues, doesn't create overlaps..
 */
void FinalExtend(MultipleAlignment *ma, double rmsd);
void CalcPartial(MultipleAlignment *ma, double rmsd);
#endif
