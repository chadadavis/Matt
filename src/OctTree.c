/* Octree.c -- Octree creation/search functions, with heap.
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

#include "OctTree.h"
#include <malloc.h>
#include <string.h>
#include <math.h>


void InitOctTreeNode(OctTreeNode *n, const Vector *minCorner, const Vector *maxCorner) {
	memset(n, 0, sizeof(OctTreeNode));
	n->corners[0] = *minCorner;
	n->corners[1] = *maxCorner;
	n->center.x = (n->corners[0].x+n->corners[1].x)/2;
	n->center.y = (n->corners[0].y+n->corners[1].y)/2;
	n->center.z = (n->corners[0].z+n->corners[1].z)/2;
}

void *GetMemory(MemoryManagementInfo *info, int size) {
	void *out;
	if (!info->numHeaps || info->heaps[info->numHeaps-1].size < size + info->heaps[info->numHeaps-1].used) {
		if (info->numHeaps % 16 == 0) {
			info->heaps = (Heap*) realloc(info->heaps, (16+info->numHeaps)*sizeof(Heap));
		}
		info->heaps[info->numHeaps].memory = (char*)malloc(1024*512);
		info->heaps[info->numHeaps].size = 1024*512;
		info->heaps[info->numHeaps++].used = 0;
	}
	out = info->heaps[info->numHeaps-1].memory+info->heaps[info->numHeaps-1].used;
	info->heaps[info->numHeaps-1].used += size;
	return out;
}

void FreeAllMemory(MemoryManagementInfo *info) {
	int i;
	for (i=0; i<info->numHeaps; i++) {
		free(info->heaps[i].memory);
	}
	free(info->heaps);
	info->heaps = 0;
	info->numHeaps = 0;
}

int GetEntries(void ** out, OctTreeNode *node, const Vector *coords, double radius) {
	int count = 0;
	Vector diff;
	subVect(&diff, coords, &node->center);
	if (lengthSquaredVect(&diff) > (node->radius+radius)*(node->radius+radius)) {
		return 0;
	}
	if (!node->numVertices) {
		if (diff.x < radius) {
			if (diff.y < radius) {
				if (node->children[0] && diff.z < radius) {
					count += GetEntries(out, node->children[0], coords, radius);
				}
				if (node->children[4] && diff.z >= -radius) {
					count += GetEntries(out+count, node->children[4], coords, radius);
				}
			}
			if (diff.y >= -radius) {
				if (node->children[2] && diff.z < radius) {
					count += GetEntries(out+count, node->children[2], coords, radius);
				}
				if (node->children[6] && diff.z >= -radius) {
					count += GetEntries(out+count, node->children[6], coords, radius);
				}
			}
		}
		if (diff.x >= -radius) {
			if (diff.y < radius) {
				if (node->children[1] && diff.z < radius) {
					count += GetEntries(out+count, node->children[1], coords, radius);
				}
				if (node->children[5] && diff.z >= -radius) {
					count += GetEntries(out+count, node->children[5], coords, radius);
				}
			}
			if (diff.y >= -radius) {
				if (node->children[3] && diff.z < radius) {
					count += GetEntries(out+count, node->children[3], coords, radius);
				}
				if (node->children[7] && diff.z >= -radius) {
					count += GetEntries(out+count, node->children[7], coords, radius);
				}
			}
		}
		return count;
	}
	else {
		int i;
		for (i=node->numVertices-1; i>=0; i--) {
			out[i] = node->vertices[i];
		}
		return node->numVertices;
	}
}


void AddEntry(OctTreeNode *node, Vector *coords, MemoryManagementInfo *info) {
	Vector diff;
	double r3;
	subVect(&diff, coords, &node->center);
	r3 = lengthSquaredVect(&diff);
	if (r3 > node->radius * node->radius) node->radius = sqrt(r3);

	if (!node->numVertices) {
		int add = (diff.x>=0) + ((diff.y>=0) << 1) + ((diff.z>=0) << 2);
		if (!node->children[add]) {
			Vector c[2];
			c[0] = node->center;
			c[1] = node->center;
			if (add & 1) {
				c[1].x = node->corners[1].x;
			}
			else {
				c[0].x = node->corners[0].x;
			}
			if (add & 2)  {
				c[1].y = node->corners[1].y;
			}
			else {
				c[0].y = node->corners[0].y;
			}

			if (add & 4)  {
				c[1].z = node->corners[1].z;
			}
			else {
				c[0].z = node->corners[0].z;
			}

			node->children[add] = (OctTreeNode*) GetMemory(info, sizeof(OctTreeNode));
			InitOctTreeNode(node->children[add], c, c+1);

			node->children[add]->numVertices = 1;
			node->children[add]->vertices[0] = coords;
			subVect(&diff, coords, &node->children[add]->center);
			node->children[add]->radius = lengthVect(&diff);
		}
		else
			AddEntry(node->children[add], coords, info);
	}
	else if (node->numVertices == MAX_OCT_TREE_VERTICES) {
		int i;
		Vector *vertices[MAX_OCT_TREE_VERTICES];
		for (i=MAX_OCT_TREE_VERTICES-1; i>=0; i--) {
			vertices[i] = node->vertices[i];
			node->vertices[i] = 0;
		}
		node->numVertices = 0;
		AddEntry(node, coords, info);
		for (i=MAX_OCT_TREE_VERTICES-1; i>=0; i--) {
			AddEntry(node, vertices[i], info);
		}
	}
	else {
		node->vertices[node->numVertices++] = coords;
	}
}
