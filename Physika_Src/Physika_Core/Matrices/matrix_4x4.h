/*
 * @file matrix_4x4.h
 * @brief 4x4 matrix.
 * @author Sheng Yang, Fei Zhu, Liyou Xu
 *
 * This file is part of Physika, a versatile physics simulation library.
 * Copyright (C) 2013- Physika Group.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0.
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#ifndef PHYSIKA_CORE_MATRICES_MATRIX_4X4_H_
#define PHYSIKA_CORE_MATRICES_MATRIX_4X4_H_

#include <glm/mat4x4.hpp>

#include "Physika_Core/Utilities/physika_assert.h"
#include "Physika_Core/Utilities/cuda_utilities.h"
#include "Physika_Core/Utilities/type_utilities.h"
#include "Physika_Core/Matrices/square_matrix.h"

namespace Physika{

template <typename Scalar, int Dim> class Vector;

/*
 * SquareMatrix<Scalar,4> are defined for C++ fundamental integer and floating-point types
 */

template <typename Scalar>
class SquareMatrix<Scalar,4>
{
public:
    CPU_GPU_FUNC_DECL SquareMatrix();
    CPU_GPU_FUNC_DECL explicit SquareMatrix(Scalar);
    CPU_GPU_FUNC_DECL SquareMatrix(Scalar x00, Scalar x01, Scalar x02, Scalar x03,
                                   Scalar x10, Scalar x11, Scalar x12, Scalar x13,
                                   Scalar x20, Scalar x21, Scalar x22, Scalar x23,
                                   Scalar x30, Scalar x31, Scalar x32, Scalar x33);
    CPU_GPU_FUNC_DECL SquareMatrix(const Vector<Scalar,4> &row1, const Vector<Scalar,4> &row2, const Vector<Scalar,4> &row3, const Vector<Scalar,4> &row4);
    
    CPU_GPU_FUNC_DECL SquareMatrix(const SquareMatrix<Scalar,4>&) = default;
    CPU_GPU_FUNC_DECL ~SquareMatrix() = default;

    CPU_GPU_FUNC_DECL  static unsigned int rows() {return 4;}
    CPU_GPU_FUNC_DECL  static unsigned int cols() {return 4;}

    CPU_GPU_FUNC_DECL Scalar& operator() (unsigned int i, unsigned int j );
    CPU_GPU_FUNC_DECL const Scalar& operator() (unsigned int i, unsigned int j) const;

    CPU_GPU_FUNC_DECL const Vector<Scalar,4> rowVector(unsigned int i) const;
    CPU_GPU_FUNC_DECL const Vector<Scalar,4> colVector(unsigned int i) const;

    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> operator+ (const SquareMatrix<Scalar,4> &) const;
    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator+= (const SquareMatrix<Scalar,4> &);
    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> operator- (const SquareMatrix<Scalar,4> &) const;
    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator-= (const SquareMatrix<Scalar,4> &);

    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator= (const SquareMatrix<Scalar,4> &) = default;

    CPU_GPU_FUNC_DECL bool operator== (const SquareMatrix<Scalar,4> &) const;
    CPU_GPU_FUNC_DECL bool operator!= (const SquareMatrix<Scalar,4> &) const;

    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> operator* (Scalar) const;
    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator*= (Scalar);

    CPU_GPU_FUNC_DECL const Vector<Scalar,4> operator* (const Vector<Scalar,4> &) const;
    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> operator* (const SquareMatrix<Scalar,4> &) const;
    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator*= (const SquareMatrix<Scalar,4> &);

    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> operator/ (Scalar) const;
    CPU_GPU_FUNC_DECL SquareMatrix<Scalar,4>& operator/= (Scalar);

    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar, 4> operator- (void) const;

    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> transpose() const;
    CPU_GPU_FUNC_DECL const SquareMatrix<Scalar,4> inverse() const;

    CPU_GPU_FUNC_DECL Scalar determinant() const;
    CPU_GPU_FUNC_DECL Scalar trace() const;
    CPU_GPU_FUNC_DECL Scalar doubleContraction(const SquareMatrix<Scalar,4> &) const;//double contraction
    CPU_GPU_FUNC_DECL Scalar frobeniusNorm() const;

    void singularValueDecomposition(SquareMatrix<Scalar,4> &left_singular_vectors,
                                    Vector<Scalar,4> &singular_values, //singular values are in descending order
                                    SquareMatrix<Scalar,4> &right_singular_vectors) const;

    void singularValueDecomposition(SquareMatrix<Scalar,4> &left_singular_vectors,
                                    SquareMatrix<Scalar,4> &singular_values_diagonal,   //singular values in descending order as a diagonal matrix
                                    SquareMatrix<Scalar,4> &right_singular_vectors) const;

    void eigenDecomposition(Vector<Scalar,4> &eigen_values_real, 
                            Vector<Scalar,4> &eigen_values_imag,
                            SquareMatrix<Scalar,4> &eigen_vectors_real, 
                            SquareMatrix<Scalar,4> &eigen_vectors_imag);

    CPU_GPU_FUNC_DECL static const SquareMatrix<Scalar,4> identityMatrix();

protected:
    glm::tmat4x4<Scalar> data_; //default: zero matrix

private:
    void compileTimeCheck()
    {
        //SquareMatrix<Scalar,Dim> is only defined for element type of integers and floating-point types
        //compile time check
        PHYSIKA_STATIC_ASSERT((is_integer<Scalar>::value || is_floating_point<Scalar>::value),
                              "SquareMatrix<Scalar,4> are only defined for integers types and floating-point types.");
    }

};

//overriding << for SquareMatrix<Scalar,4>
template <typename Scalar>
inline std::ostream& operator<< (std::ostream &s, const SquareMatrix<Scalar,4> &mat)
{
    if((is_same<Scalar,unsigned char>::value)||(is_same<Scalar,signed char>::value))
    {
        s<<"[";
        s<<static_cast<int>(mat(0,0))<<", "<<static_cast<int>(mat(0,1))<<", "<<static_cast<int>(mat(0,2))<<", "<<static_cast<int>(mat(0,3))<<"; ";
        s<<static_cast<int>(mat(1,0))<<", "<<static_cast<int>(mat(1,1))<<", "<<static_cast<int>(mat(1,2))<<", "<<static_cast<int>(mat(1,3))<<"; ";
        s<<static_cast<int>(mat(2,0))<<", "<<static_cast<int>(mat(2,1))<<", "<<static_cast<int>(mat(2,2))<<", "<<static_cast<int>(mat(2,3))<<"; ";
        s<<static_cast<int>(mat(3,0))<<", "<<static_cast<int>(mat(3,1))<<", "<<static_cast<int>(mat(3,2))<<", "<<static_cast<int>(mat(3,3))<<"]";
    }
    else
    {
        s<<"[";
        s<<mat(0,0)<<", "<<mat(0,1)<<", "<<mat(0,2)<<", "<<mat(0,3)<<"; ";
        s<<mat(1,0)<<", "<<mat(1,1)<<", "<<mat(1,2)<<", "<<mat(1,3)<<"; ";
        s<<mat(2,0)<<", "<<mat(2,1)<<", "<<mat(2,2)<<", "<<mat(2,3)<<"; ";
        s<<mat(3,0)<<", "<<mat(3,1)<<", "<<mat(3,2)<<", "<<mat(3,3)<<"]";
    }
    return s;
}

//make * operator commutative
template <typename S, typename T>
CPU_GPU_FUNC_DECL  const SquareMatrix<T,4> operator* (S scale, const SquareMatrix<T,4> &mat)
{
    return mat*scale;
}

//convenient typedefs
typedef SquareMatrix<float,4> Matrix4f;
typedef SquareMatrix<double,4> Matrix4d;
typedef SquareMatrix<int,4> Matrix4i;

}  //end of namespace Physika

#endif //PHYSIKA_CORE_MATRICES_MATRIX_4X4_H_
