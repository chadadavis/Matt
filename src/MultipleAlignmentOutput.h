/*
 *  Contains all functions to output an alignment to a file.
 */

#ifndef MULTIPLE_ALIGNMENT_OUTPUT_H
#define MULTIPLE_ALIGNMENT_OUTPUT_H

#include "MultipleAlignment.h"
#include <stdio.h>

/* Order is only needed to get the right base structure for the bent function.
 * Order has ids, not indices.
 */
int OutputAlignmentBentPDB(MultipleAlignment *ma, const char *fileName, int *order);
int OutputAlignmentPDB(MultipleAlignment *ma, const char *filename, int *order);
int OutputAlignmentUnbentPDB(MultipleAlignment *ma, const char *filename);

int OutputAlignmentRasmol(MultipleAlignment *ma, const char *filename);
int OutputAlignmentSpdbv(MultipleAlignment *ma, const char *filename);

int OutputAlignmentText(MultipleAlignment *ma, const char *filename, int root);
int OutputAlignmentFasta(MultipleAlignment *ma, const char *filename);
int OutputAlignmentMsf(MultipleAlignment *ma, char *filename);
int OutputAlignmentPir(MultipleAlignment *ma, const char *filename);

/* Not strictly an output function, but only used for getting an alignment to output. */
void ReorderAlignment(MultipleAlignment *ma, int *order);
#endif
