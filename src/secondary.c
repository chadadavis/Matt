/* Protein.h -- Functions for dealing with and detecting alpha helices and
 * beta sheets.
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

#include "secondary.h"

int FindFirstMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair) {
	if (betaPair->start1 <= pos && pos < betaPair->start1+betaPair->length) {
		residuePair->res1 = pos;
		residuePair->res2 = betaPair->start2+(pos-betaPair->start1)*betaPair->direction;
		residuePair->dir = betaPair->direction;
		return 1;
	}
	return 0;
}

int FindSecondMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair) {
	if (betaPair->direction == 1) {
		if (betaPair->start2 <= pos && pos < betaPair->start2+betaPair->length) {
			residuePair->res1 = pos;
			residuePair->res2 = betaPair->start1+(pos-betaPair->start2);
			residuePair->dir = betaPair->direction;
			return 1;
		}
	}
	else {
		if (betaPair->start2 >= pos && pos > betaPair->start2-betaPair->length) {
			residuePair->res1 = pos;
			residuePair->res2 = betaPair->start1-(pos-betaPair->start2);
			residuePair->dir = betaPair->direction;
			return 1;
		}
	}
	return 0;
}

void FlipPair(BetaPair *betaPair) {
	if (betaPair->direction == 1) {
		int temp = betaPair->start1;
		betaPair->start1 = betaPair->start2;
		betaPair->start2 = temp;
	}
	else {
		int temp = betaPair->start1 + betaPair->length - 1;
		betaPair->start1 = betaPair->start2 - betaPair->length + 1;
		betaPair->start2 = temp;
	}
}

int FindMatch(BetaPair *betaPair, int pos, BetaResiduePair *residuePair) {
	if (FindFirstMatch(betaPair, pos, residuePair)) return 1;
	return FindSecondMatch(betaPair, pos, residuePair);
}

void FlipSheet(BetaSheet *b) {
	int i;
	for (i = b->numStrands/2-1; i>=0; i--) {
		BetaSheetPair temp = b->strands[i];
		b->strands[i] = b->strands[b->numStrands-1-i];
		b->strands[b->numStrands-1-i] = temp;
	}
	for (i = b->numStrands-1; i>=0; i--) {
		if (b->strands[i].direction == 1) {
			int temp = b->strands[i].res1;
			b->strands[i].res1 = b->strands[i].res2;
			b->strands[i].res2 = temp;
		}
		else {
			int temp = b->strands[i].res1 + b->strands[i].length - 1;
			b->strands[i].res1 = b->strands[i].res2 - b->strands[i].length + 1;
			b->strands[i].res2 = temp;
		}
	}
}
