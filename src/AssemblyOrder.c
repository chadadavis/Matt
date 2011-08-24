/* AssemblyOrder.h -- Tree creation functions.
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

#include "AssemblyOrder.h"
#include <malloc.h>

AssemblyOrder *MakeAssemblyOrder(int id) {
	AssemblyOrder *out = (AssemblyOrder*) malloc(sizeof(AssemblyOrder) + sizeof(AssemblyNode));
	out->nodes = 1;
	out->root = (AssemblyNode*) (out+1);
	out->root->right = out->root->left = 0;
	out->root->id = id;
	return out;
};

void CopyAssemblyNode(AssemblyNode *dest, AssemblyNode *source, AssemblyNode** last) {
	dest->id = source->id;
	dest->left = source->left;
	if (source->left) {
		dest->left = ++last[0];
		CopyAssemblyNode(dest->left, source->left, last);
	}
	dest->right = source->right;
	if (source->right) {
		dest->right = ++last[0];
		CopyAssemblyNode(dest->right, source->right, last);
	}
};

AssemblyOrder *DuplicateAssemblyOrder(AssemblyOrder *order) {
	AssemblyOrder *out = (AssemblyOrder*) malloc(sizeof(AssemblyOrder) + sizeof(AssemblyNode)*order->nodes);
	AssemblyNode *last = out->root = (AssemblyNode*) (out+1);
	out->nodes = order->nodes;
	CopyAssemblyNode(out->root, order->root, &last);
	return out;
};

AssemblyOrder *CombineAssemblyOrders(AssemblyOrder *left, AssemblyOrder *right) {
	int nodes = left->nodes + right->nodes + 1;
	AssemblyOrder *out = (AssemblyOrder*) malloc(sizeof(AssemblyOrder) + sizeof(AssemblyNode)*nodes);
	AssemblyNode *last = out->root = (AssemblyNode*) (out+1);
	out->nodes = nodes;
	out->root->left = out->root->right = 0;
	out->root->id = -1;
	out->root->left = ++last;
	CopyAssemblyNode(out->root->left, left->root, &last);
	out->root->right = ++last;
	CopyAssemblyNode(out->root->right, right->root, &last);
	return out;
};
