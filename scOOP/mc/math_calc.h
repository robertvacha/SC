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
    arg+=6755399441055744.0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
    return static_cast<double>(reinterpret_cast<int&>(arg) );
#pragma GCC diagnostic pop
    /*if (arg < 0) {
        return (double)( (long)(arg-0.5) );
    } else {
        return (double)( (long)(arg+0.5) );
    }*/
}


/**
 * @brief vec cross_product
 * @param A
 * @param B
 * @return
 */
inline Vector vecCrossProduct(const Vector& A, const Vector& B) {
    return Vector(A.y*B.z - A.z*B.y, -A.x*B.z + A.z*B.x, A.x*B.y - A.y*B.x);
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
