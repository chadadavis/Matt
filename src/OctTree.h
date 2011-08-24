/* Octree.h -- Octree creation/search prototypes, with heap.
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

#ifndef OCT_TREE_H
#define OCT_TREE_H

#include "Vector.h"

#define MAX_OCT_TREE_VERTICES 6

struct OCT_TREE_NODE {
	Vector center;
	Vector corners[2];

	union {
		struct OCT_TREE_NODE *children[8];
		Vector *vertices[MAX_OCT_TREE_VERTICES];
	};
	double radius;
	int numVertices;
};
typedef struct OCT_TREE_NODE OctTreeNode;

struct HEAP {
	int size;
	int used;
	char *memory;
};
typedef struct HEAP Heap;

struct MEMORY_MANAGEMENT_INFO {
	int numHeaps;
	Heap *heaps;
};
typedef struct MEMORY_MANAGEMENT_INFO MemoryManagementInfo;

void *GetMemory(MemoryManagementInfo *info, int size);
void FreeAllMemory(MemoryManagementInfo *info);

void InitOctTreeNode(OctTreeNode *n, const Vector *minCorner, const Vector *maxCorner);
void AddEntry(OctTreeNode *node, Vector *coords, MemoryManagementInfo *info);
int GetEntries(void ** out, OctTreeNode *node, const Vector *coords, double radius);

#endif
