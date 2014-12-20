/** @file math_calc.h*/

#ifndef MATH_CALC2_H
#define MATH_CALC2_H

#include "../structures/Vector.h"


/**
 * @brief Returns the nearest integer to its argument as a double precision number. e.g.
   anint(-0.49) = 0.0 and anint(-0.51) = -1.0. Equivalent to the Fortran intrinsic ANINT.
 * @param arg
 * @return
 */
inline double anInt(double arg)  {
    if (arg < 0) {
        return (double)( (long)(arg-0.5) );
    } else {
        return (double)( (long)(arg+0.5) );
    }
}

/**
 * @brief Returns the vector pointing from the centre of mass of particle 2 to the
   centre of mass of the closest image of particle 1.
 * @param r1
 * @param r2
 * @param box
 * @return
 */
inline Vector image(Vector* r1, Vector* r2, Vector* box) {
    double x = r1->x - r2->x,y = r1->y - r2->y,z = r1->z - r2->z;
    return Vector( box->x * (x - anInt(x)),
                  box->y * (y - anInt(y)),
                  box->z * (z - anInt(z)) );
}


/**
 * @brief vec cross_product
 * @param A
 * @param B
 * @return
 */
inline Vector vecCrossProduct(Vector* A, Vector* B) {
    return Vector(A->y*B->z - A->z*B->y, -A->x*B->z + A->z*B->x, A->x*B->y - A->y*B->x);
}

/**
 * @brief ortogonalise
 * @param A
 * @param B
 */
inline void ortogonalise(Vector *A, Vector *B) {
    double dp(A->dot(*B));    A->x -= dp*B->x;    A->y -= dp*B->y;    A->z -= dp*B->z;
}


#endif // MATH_CALC2_H
