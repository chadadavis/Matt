/* AssemblyOrder.h -- Prototypes for creating a tree from two trees and a new root node.
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

#ifndef ASSEMBLY_ORDER_H
#define ASSEMBLY_ORDER_H

struct ASSEMBLY_NODE {
	int id;
	struct ASSEMBLY_NODE *right, *left;
};
typedef struct ASSEMBLY_NODE AssemblyNode;

struct ASSEMBLY_ORDER {
	int nodes;
	AssemblyNode *root;
};
typedef struct ASSEMBLY_ORDER AssemblyOrder;


AssemblyOrder *MakeAssemblyOrder(int id);
void CopyAssemblyNode(AssemblyNode *dest, AssemblyNode *source, AssemblyNode** last);
AssemblyOrder *DuplicateAssemblyOrder(AssemblyOrder *order);
AssemblyOrder *CombineAssemblyOrders(AssemblyOrder *left, AssemblyOrder *right);

#endif
