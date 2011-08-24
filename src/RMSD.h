/* RMSD.h -- RMSD calculation prototypes.
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

#include "Vector.h"

double CalcRMSDAndRotation(Matrix *rotation, Matrix *A, double E0, int len);
/* ATA is written to for calling calc rotation, if needed. */
double CalcRMSD(Matrix *ATA, Matrix *A, double E0, int len);
double CalcRotation(Matrix *rotation, Matrix *A, Matrix *ATA, double E0, int len);
