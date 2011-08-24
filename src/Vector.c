/* Vector.c -- Vector, Matrix, and Quaternion computation functions.
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
#include <math.h>

double detMat(Matrix *m) {
	return
		m->M[0][0] * m->M[1][1] * m->M[2][2] -
		m->M[0][0] * m->M[1][2] * m->M[2][1] +
		m->M[1][0] * m->M[2][1] * m->M[0][2] -
		m->M[2][2] * m->M[1][2] * m->M[2][1] +
		m->M[2][0] * m->M[1][2] * m->M[0][1] -
		m->M[1][1] * m->M[0][2] * m->M[2][0];
}

Vector *subVect(Vector *v1, const Vector *v2, const Vector *v3) {
	v1->x = v2->x - v3->x;
	v1->y = v2->y - v3->y;
	v1->z = v2->z - v3->z;
	return v1;
}

Vector *addVect(Vector *v1, const Vector *v2, const Vector *v3) {
	v1->x = v2->x + v3->x;
	v1->y = v2->y + v3->y;
	v1->z = v2->z + v3->z;
	return v1;
}

Vector *mulVect(Vector *v1, const Vector *v2, double f) {
	v1->x = v2->x * f;
	v1->y = v2->y * f;
	v1->z = v2->z * f;
	return v1;
}

Vector *divVect(Vector *v1, const Vector *v2, double f) {
	v1->x = v2->x / f;
	v1->y = v2->y / f;
	v1->z = v2->z / f;
	return v1;
}

double dotVect(const Vector *v1, const Vector *v2) {
	return v2->x*v1->x + v2->y*v1->y + v2->z*v1->z;
}

Vector *crossVect(Vector *v1, const Vector *v2, const Vector *v3) {
	v1->x = v2->y * v3->z - v2->z * v3->y;
	v1->y = v2->z * v3->x - v2->x * v3->z,
	v1->z = v2->x * v3->y - v2->y * v3->x;
	return v1;
}

double lengthSquaredVect(const Vector *v) {
	return dotVect(v, v);
}

double lengthVect(const Vector *v) {
	return sqrt(dotVect(v,v));
}

Vector *normalizeVect(Vector *v) {
	return divVect(v, v, lengthVect(v));
}

double distVect(Vector *v1, Vector *v2) {
	double dx = v1->x - v2->x;
	double dy = v1->y - v2->y;
	double dz = v1->z - v2->z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

double distSquaredVect(Vector *v1, Vector *v2) {
	double dx = v1->x - v2->x;
	double dy = v1->y - v2->y;
	double dz = v1->z - v2->z;
	return dx*dx+dy*dy+dz*dz;
}


const Matrix identity = {1,0,0,0,0,1,0,0,0,0,1,0};

void mulMat(Matrix *out, const Matrix *in1, const Matrix *in2) {
	int row, col, i;
	for (row = 0; row<3; row++) {
		out->M[row][0] = 0;
		out->M[row][1] = 0;
		out->M[row][2] = 0;
		out->M[row][3] = in1->M[row][3];
		for (col = 0; col<4; col++) {
			for (i=0; i<3; i++) {
				out->M[row][col] += in1->M[row][i] * in2->M[i][col];
			}
		}
	}
}

void mulMat3x3(Matrix *out, const Matrix *in1, const Matrix *in2) {
	int row, col;
	for (row = 0; row<3; row++) {
		for (col = 0; col<3; col++) {
			out->M[row][col] = in1->M[row][0] * in2->M[0][col] +
				in1->M[row][1] * in2->M[1][col] +
				in1->M[row][2] * in2->M[2][col];
		}
	}
}

void transformVect(Vector *out, const Matrix *m, const Vector *in) {
	out->x = m->M[0][0] * in->x + m->M[0][1] * in->y + m->M[0][2] * in->z + m->M[0][3];
	out->y = m->M[1][0] * in->x + m->M[1][1] * in->y + m->M[1][2] * in->z + m->M[1][3];
	out->z = m->M[2][0] * in->x + m->M[2][1] * in->y + m->M[2][2] * in->z + m->M[2][3];
}

void invertMat(Matrix *out, const Matrix *in) {
	out->M[0][0] = in->M[1][1]*in->M[2][2] - in->M[2][1]*in->M[1][2];
	out->M[1][0] = in->M[2][0]*in->M[1][2] - in->M[1][0]*in->M[2][2];
	out->M[2][0] = in->M[1][0]*in->M[2][1] - in->M[2][0]*in->M[1][1];

	out->M[0][1] = in->M[2][1]*in->M[0][2] - in->M[0][1]*in->M[2][2];
	out->M[1][1] = in->M[0][0]*in->M[2][2] - in->M[2][0]*in->M[0][2];
	out->M[2][1] = in->M[2][0]*in->M[0][1] - in->M[0][0]*in->M[2][1];

	out->M[0][2] = in->M[0][1]*in->M[1][2] - in->M[1][1]*in->M[0][2];
	out->M[1][2] = in->M[1][0]*in->M[0][2] - in->M[0][0]*in->M[1][2];
	out->M[2][2] = in->M[0][0]*in->M[1][1] - in->M[1][0]*in->M[0][1];

	out->M[0][3] = in->M[0][1]*(in->M[1][3]*in->M[2][2] - in->M[1][2]*in->M[2][3]) + in->M[0][2]*(in->M[1][1]*in->M[2][3] - in->M[1][3]*in->M[2][1]) + in->M[0][3]*(in->M[1][2]*in->M[2][1] - in->M[1][1]*in->M[2][2]);
	out->M[1][3] = in->M[0][0]*(in->M[1][2]*in->M[2][3] - in->M[1][3]*in->M[2][2]) + in->M[0][2]*(in->M[1][3]*in->M[2][0] - in->M[1][0]*in->M[2][3]) + in->M[0][3]*(in->M[1][0]*in->M[2][2] - in->M[1][2]*in->M[2][0]);
	out->M[2][3] = in->M[0][1]*(in->M[1][0]*in->M[2][3] - in->M[1][3]*in->M[2][0]) + in->M[0][0]*(in->M[1][3]*in->M[2][1] - in->M[1][1]*in->M[2][3]) + in->M[0][3]*(in->M[1][1]*in->M[2][0] - in->M[1][0]*in->M[2][1]);

}

void matrixToQuaternion(Quaternion *out, const Matrix *in) {
	double trace = in->M[0][0] + in->M[1][1] + in->M[2][2] + 1;
	double S;
	if (trace > 0.0000001) {
		S = sqrt(trace);
		out->w = 0.5 * S;
		S = 0.5 / S;
		out->v.x = (in->M[2][1] - in->M[1][2]) * S;
		out->v.y = (in->M[0][2] - in->M[2][0]) * S;
		out->v.z = (in->M[1][0] - in->M[0][1]) * S;
	}
	else if (in->M[0][0] > in->M[1][1] && in->M[0][0] > in->M[2][2]) {
		S = 2*sqrt(1 + in->M[0][0] - in->M[1][1] - in->M[2][2]);
		out->v.x = 0.25 * S;
		out->v.y = (in->M[1][0] + in->M[0][1]) / S;
		out->v.z = (in->M[2][0] + in->M[0][2]) / S;
		out->w = (in->M[1][2] - in->M[2][1]) / S;
	}
	else if (in->M[1][1] > in->M[2][2]) {
		S = 2*sqrt(1 - in->M[0][0] + in->M[1][1] - in->M[2][2]);
		out->v.x = (in->M[1][0] + in->M[0][1]) / S;
		out->v.y = 0.25 * S;
		out->v.z = (in->M[2][1] + in->M[1][2]) / S;
		out->w = (in->M[2][0] - in->M[0][2]) / S;
	}
	else {
		S = 2*sqrt(1 - in->M[0][0] - in->M[1][1] + in->M[2][2]);
		out->v.x = (in->M[2][0] + in->M[0][2]) / S;
		out->v.y = (in->M[2][1] + in->M[1][2]) / S;
		out->v.z = 0.25 * S;
		out->w = (in->M[0][1] - in->M[1][0]) / S;
	}
}

double invHalfCosQuat(const Quaternion *q1, const Quaternion *q2) {
	return fabs(q1->w*q2->w + dotVect(&q1->v, &q2->v));
}
