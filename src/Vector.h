/* Vector.h -- Vector, Matrix, and Quaternion computation prototypes.
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

#ifndef VECTOR_H
#define VECTOR_H


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Simple 3 element Vector class */
struct VECTOR {
	union {
		struct {
			double x,y,z;
		};
		double v[3];
	};
};
typedef struct VECTOR Vector;

Vector *subVect(Vector *v1, const Vector *v2, const Vector *v3);
Vector *addVect(Vector *v1, const Vector *v2, const Vector *v3);
Vector *mulVect(Vector *v1, const Vector *v2, double f);
Vector *divVect(Vector *v1, const Vector *v2, double f);
double dotVect(const Vector *v1, const Vector *v2);
double lengthSquaredVect(const Vector *v);
double lengthVect(const Vector *v);
Vector *normalizeVect(Vector *v);
double distVect(Vector *v1, Vector *v2);
double distSquaredVect(Vector *v1, Vector *v2);
/* v1 cannot be v2 or v3, unlike above functions.
 */
Vector *crossVect(Vector *v1, const Vector *v2, const Vector *v3);

struct MATRIX {
	/* First value is column (x), second is row (y).
	 * Last row is always to be {0, 0, 0, 1}, so
	 * to save memory, not stored.
	 */
	double M[3][4];
};
typedef struct MATRIX Matrix;

/* out must be distinct from input pointers, unlike in purely vector functions.
 */
void transformVect(Vector *out, const Matrix *m, const Vector *in);



void mulMat(Matrix *out, const Matrix *in1, const Matrix *in2);
void mulMat3x3(Matrix *out, const Matrix *in1, const Matrix *in2);
void invertMat(Matrix *out, const Matrix *in);

/* Identity matrix.  Nice to have handy.
 */
extern const Matrix identity;

struct QUATERNION {
	Vector v;
	double w;
};
typedef struct QUATERNION Quaternion;

/* Creates a rotation matrix that represents the
 * same rotation.
 */
void quatToMat(Matrix *out, const Quaternion *in);
void matrixToQuaternion(Quaternion *out, const Matrix *in);
double invHalfCosQuat(const Quaternion *q1, const Quaternion *q2);

double detMat(Matrix *m);

void printMat(const Matrix* m);

#endif
