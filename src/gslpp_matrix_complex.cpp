/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#ifndef __GSL_BLAS_H__
#include <gsl/gsl_blas.h>
#endif
#ifndef __GSL_LINALG_H__
#include <gsl/gsl_linalg.h>
#endif
#ifndef __GSL_EIGEN_H__
#include <gsl/gsl_eigen.h>
#endif
#ifndef GSLPP_MATRIX_COMPLEX_H
#include "gslpp_matrix_complex.h"
#endif

namespace gslpp {

    /** Constructor  */
    matrix<complex>::matrix(const size_t& size_i, const size_t& size_j, const complex& z)
    {
        _matrix = gsl_matrix_complex_alloc(size_i, size_j);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    }

    matrix<complex>::matrix(const size_t& size_i, const complex& z)
    {
        _matrix = gsl_matrix_complex_alloc(size_i, size_i);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    }

    matrix<complex>::matrix(const size_t& size_i, const size_t& size_j, const double& a)
    {
        complex z(a);
        _matrix = gsl_matrix_complex_alloc(size_i, size_j);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    }

    matrix<complex>::matrix(const size_t& size_i, const double& a)
    {
        complex z(a);
        _matrix = gsl_matrix_complex_alloc(size_i, size_i);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
    }

    /** Copy constructor  */
    matrix<complex>::matrix(const matrix<complex>& m)
    {
        _matrix = gsl_matrix_complex_alloc(m.size_i(), m.size_j());
        gsl_matrix_complex_memcpy(_matrix, m.as_gsl_type_ptr());
    }

    matrix<complex>::matrix(const matrix<double>& m)
    {
        size_t i, j, size_i, size_j;
        size_i = m.size_i();
        size_j = m.size_j();
        _matrix = gsl_matrix_complex_alloc(size_i, size_j);
        for (i = 0; i < size_i; ++i)
            for (j = 0; j < size_j; ++j)
                gsl_matrix_complex_set(_matrix, i, j, m(i, j) + 0. * complex::i());
    }

    matrix<complex>::matrix(const vector<complex>& v)
    {
        size_t i, size = v.size();
        complex z(0, 0, false);
        _matrix = gsl_matrix_complex_alloc(size, size);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
        for (i = 0; i < size; ++i)
            gsl_matrix_complex_set(_matrix, i, i, v(i));
    }

    matrix<complex>::matrix(const vector<double>& v)
    {
        size_t i, size = v.size();
        complex z(0, 0, false);
        _matrix = gsl_matrix_complex_alloc(size, size);
        gsl_matrix_complex_set_all(_matrix, z.as_gsl_type());
        for (i = 0; i < size; ++i)
            gsl_matrix_complex_set(_matrix, i, i, v(i) + 0. * complex::i());
    }

    matrix<complex>::matrix(const gsl_matrix_complex& m)
    {
        _matrix = gsl_matrix_complex_alloc(m.size1, m.size2);
        gsl_matrix_complex_memcpy(_matrix, &m);
    }

    matrix<complex>::matrix(const gsl_matrix_complex* m)
    {
        _matrix = gsl_matrix_complex_alloc(m->size1, m->size2);
        gsl_matrix_complex_memcpy(_matrix, m);
    }

    /** Destructor */
    matrix<complex>::~matrix()
    {
        gsl_matrix_complex_free(_matrix);
    }

    /** Get i-th element */
    const complex matrix<complex>::operator()(const size_t& i, const size_t& j) const
    {
        gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
        return complex(x);
    }

    /** Set i-th element */
    //  complex& matrix<complex>::operator()(const size_t& i, const size_t& j)
    //  {
    //    gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
    //    return complex(x);
    //  }

    /** Assign */
  void matrix<complex>::reset()
  {
      gsl_matrix_complex_set_zero(_matrix);
  }
  
    matrix<complex>& matrix<complex>::operator=(const matrix<complex>& m)
    {
        if (size_i() == m.size_i() && size_j() == m.size_j())
            gsl_matrix_complex_memcpy(_matrix, m.as_gsl_type_ptr());
        else
        {
            std::cout << "\n ** Wrong size assign in matrix<complex> operator =" << std::endl;
            exit(EXIT_FAILURE);
        }
        return *this;
    }

    /** Assign element */
    void matrix<complex>::assign(const size_t& i, const size_t& j, const complex& z)
    {
        gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
        *x = z.as_gsl_type();
    }

    void matrix<complex>::assign(const size_t& i, const size_t& j, const double& a)
    {
        gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
        *x = complex(a).as_gsl_type();
    }

    void matrix<complex>::assignre(const size_t& i, const size_t& j, const double& a)
    {
        gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
        GSL_SET_REAL(x, a);
    }

    void matrix<complex>::assignim(const size_t& i, const size_t& j, const double& a)
    {
        gsl_complex *x = gsl_matrix_complex_ptr(_matrix, i, j);
        GSL_SET_IMAG(x, a);
    }

    /** Assign submatrix */
    void matrix<complex>::assign(const size_t& i, const size_t& j, const matrix<complex>& z)
    {
        size_t ki, kj;
        gsl_complex *x;
        if (i + z.size_i() <= size_i() && j + z.size_j() <= size_j())
            for (ki = i; ki < i + z.size_i(); ki++)
                for (kj = j; kj < j + z.size_j(); kj++)
                {
                    x = gsl_matrix_complex_ptr(_matrix, ki, kj);
                    *x = z(ki - i, kj - j).as_gsl_type();
                } else
        {
            std::cout << "\n ** Wrong size assign in matrix<complex> assign submatrix" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void matrix<complex>::assign(const size_t& i, const size_t& j, const matrix<double>& a)
    {
        matrix<complex> z(a);
        assign(i, j, z);
    }

    /** Get matrix size */
    size_t matrix<complex>::size_i() const
    {
        return _matrix->size1;
    }

    size_t matrix<complex>::size_j() const
    {
        return _matrix->size2;
    }

    size_t matrix<complex>::size() const
    {
        if (_matrix->size1 == _matrix->size2)
            return _matrix->size2;
        else
        {
            std::cout << "\n ** Non square matrix" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    /** Identity matrix */
    matrix<complex> matrix<complex>::Id(size_t size)
    {
        return matrix<double>::Id(size) * (1. + 0. * complex::i());
    }

    /** Transpose matrix */
    matrix<complex> matrix<complex>::transpose() const
    {
        matrix<complex> m1(size_j(), size_i(), 0.);
        if (gsl_matrix_complex_transpose_memcpy(m1.as_gsl_type_ptr(), _matrix))
        {
            std::cout << "\n ** Error in matrix<complex> transpose" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Hermitean conjugate matrix */
    matrix<complex> matrix<complex>::hconjugate() const
    {
        matrix<complex> m1(size_j(), size_i(), 0.);
        for (size_t i = 0; i < size_i(); ++i)
            for (size_t j = 0; j < size_j(); ++j)
                m1.assign(j, i, (*this)(i, j).conjugate());
        return m1;
    }

    /** Inverse matrix */
    matrix<complex> matrix<complex>::inverse() const
    {
        matrix<complex> m1(_matrix);
        matrix<complex> m2(_matrix);
        int signum;
        gsl_permutation *p;

        if (size_j() != size_i())
        {
            std::cout << "\n ** Size mismatch in matrix<complex> inverse" << std::endl;
            exit(EXIT_FAILURE);
        }

        if ((p = gsl_permutation_alloc(size_i())) == NULL)
        {
            std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (gsl_linalg_complex_LU_decomp(m1.as_gsl_type_ptr(), p, &signum))
        {
            std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
            gsl_permutation_free(p);
            exit(EXIT_FAILURE);
        }

        if (gsl_linalg_complex_LU_invert(m1.as_gsl_type_ptr(), p,
                m2.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> inverse" << std::endl;
            gsl_permutation_free(p);
            exit(EXIT_FAILURE);
        }
        gsl_permutation_free(p);
        return m2;
    }
    
    /** Check if matrix is singular (This is to avoid error checking by luc.c)**/

    bool matrix<complex>::isSingular()
    {
        matrix<complex> m1(_matrix);
        int signum;
        gsl_permutation *p;

        if (size_j() != size_i())
        {
            std::cout << "\n ** Size mismatch in bool isSingular" << std::endl;
            exit(EXIT_FAILURE);
        }

        if ((p = gsl_permutation_alloc(size_i())) == NULL)
        {
            std::cout << "\n ** Error in bool isSingular" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (gsl_linalg_complex_LU_decomp(m1.as_gsl_type_ptr(), p, &signum))
        {
            std::cout << "\n ** Error in bool isSingular" << std::endl;
            gsl_permutation_free(p);
            exit(EXIT_FAILURE);
        }
        
        size_t i, n = m1.size_i();

        for (i = 0; i < n; i++) {
            gsl_complex u = gsl_matrix_complex_get(m1.as_gsl_type_ptr(), i, i);
            if (GSL_REAL(u) == 0 && GSL_IMAG(u) == 0) return 1;
        }
        return 0;
    }

    /** Get matrix of real parts */
    matrix<double> matrix<complex>::real() const
    {
        matrix<double> res(size_i(), size_j());
        for (size_t i = 0; i < size_i(); i++)
        {
            gsl_vector_complex_view tmp = gsl_matrix_complex_row(_matrix, i);
            gsl_vector_view tmpd = gsl_vector_complex_real(&tmp.vector);
            gsl_matrix_set_row(res.as_gsl_type_ptr(), i, &tmpd.vector);
        }
        return res;
    }

    /** Get matrix of imaginary parts */
    matrix<double> matrix<complex>::imag() const
    {
        matrix<double> res(size_i(), size_j());
        for (size_t i = 0; i < size_i(); i++)
        {
            gsl_vector_complex_view tmp = gsl_matrix_complex_row(_matrix, i);
            gsl_vector_view tmpd = gsl_vector_complex_imag(&tmp.vector);
            gsl_matrix_set_row(res.as_gsl_type_ptr(), i, &tmpd.vector);
        }
        return res;
    }

    void matrix<complex>::eigensystem(matrix<complex> &U, vector<double> &S) const
    {
        matrix<complex> m(*this);

        gsl_eigen_hermv_workspace *ws;

        ws = gsl_eigen_hermv_alloc(size_i());

        gsl_eigen_hermv(m.as_gsl_type_ptr(), S.as_gsl_type_ptr(), U.as_gsl_type_ptr(), ws);

        gsl_eigen_hermv_sort(S.as_gsl_type_ptr(), U.as_gsl_type_ptr(),
                GSL_EIGEN_SORT_VAL_ASC);

        gsl_eigen_hermv_free(ws);
    }

    namespace {
        using WideComplex = std::complex<long double>;
        using Column = std::vector<WideComplex>;

        struct JacobiColumn {
            explicit JacobiColumn(size_t size) : left(size), right(size), roundoff(size) {}

            Column left;
            Column right;
            std::vector<long double> roundoff;
            long double norm = 0.;
        };

        long double columnNorm(const Column& column)
        {
            long double norm = 0.;
            for (const auto& entry : column)
                norm = std::hypot(norm, std::abs(entry));
            return norm;
        }

        WideComplex columnInnerProduct(const Column& x, const Column& y,
                long double scaleX = 1., long double scaleY = 1.)
        {
            WideComplex product;
            for (size_t row = 0; row < x.size(); ++row)
                product += std::conj(x[row] / scaleX) * (y[row] / scaleY);
            return product;
        }

        void rotateEntries(WideComplex& x, WideComplex& y, long double cosine,
                const WideComplex& sine, const WideComplex& conjugateSine)
        {
            const auto oldX = x;
            x = cosine * oldX - conjugateSine * y;
            y = sine * oldX + cosine * y;
        }

        bool jacobiRotate(JacobiColumn& x, JacobiColumn& y, long double tolerance)
        {
            if (x.norm == 0. || y.norm == 0.) return false;

            const auto overlap = columnInnerProduct(x.left, y.left, x.norm, y.norm);
            const long double magnitude = std::abs(overlap);
            if (magnitude <= tolerance) return false;

            // Scale the Hermitian 2x2 rotation by the larger column norm squared.
            const long double ratio = std::min(x.norm, y.norm) / std::max(x.norm, y.norm);
            const long double gap = (1. - ratio) * (1. + ratio);
            const long double coupling = 2. * magnitude * ratio;
            const long double tangent = std::copysign(coupling / (gap + std::hypot(gap, coupling)),
                    y.norm - x.norm);
            if (tangent == 0.)
                throw std::runtime_error("Complex SVD rotation is outside the representable range");

            const long double cosine = 1. / std::hypot(1.L, tangent);
            const long double realSine = cosine * tangent;
            const long double sineMagnitude = std::abs(realSine);
            const WideComplex sine = realSine * (overlap / magnitude);
            const auto conjugateSine = std::conj(sine);
            const long double relativeRoundoff = 8. * std::numeric_limits<long double>::epsilon();
            for (size_t row = 0; row < x.left.size(); ++row) {
                // Estimate cancellation error for each entry from its contributing
                // entries. A small independent component can carry a singular value.
                const long double errorX = x.roundoff[row];
                const long double errorY = y.roundoff[row];
                const long double magnitudeX = std::abs(x.left[row]);
                const long double magnitudeY = std::abs(y.left[row]);
                x.roundoff[row] = std::hypot(cosine * errorX, sineMagnitude * errorY)
                        + relativeRoundoff * (cosine * magnitudeX + sineMagnitude * magnitudeY);
                y.roundoff[row] = std::hypot(sineMagnitude * errorX, cosine * errorY)
                        + relativeRoundoff * (sineMagnitude * magnitudeX + cosine * magnitudeY);

                rotateEntries(x.left[row], y.left[row], cosine, sine, conjugateSine);
                if (std::abs(x.left[row]) <= x.roundoff[row]) x.left[row] = WideComplex();
                if (std::abs(y.left[row]) <= y.roundoff[row]) y.left[row] = WideComplex();
                rotateEntries(x.right[row], y.right[row], cosine, sine, conjugateSine);
            }
            x.norm = columnNorm(x.left);
            y.norm = columnNorm(y.left);
            return true;
        }

        void completeLeftBasis(std::vector<JacobiColumn>& columns)
        {
            const size_t size = columns.size();
            // Choose projected coordinate vectors with enough norm to normalize safely.
            const long double minResidual = 1. / std::sqrt(2.L * size);
            size_t nextAxis = 0;
            Column candidate(size);
            for (auto& column : columns) {
                if (column.norm != 0.) continue;

                long double norm = 0.;
                while (nextAxis < size && norm < minResidual) {
                    std::fill(candidate.begin(), candidate.end(), WideComplex());
                    candidate[nextAxis++] = 1.;
                    for (int pass = 0; pass < 2; ++pass) {
                        for (const auto& basis : columns) {
                            const auto projection = columnInnerProduct(basis.left, candidate);
                            for (size_t row = 0; row < size; ++row)
                                candidate[row] -= projection * basis.left[row];
                        }
                    }
                    norm = columnNorm(candidate);
                }
                if (norm < minResidual)
                    throw std::runtime_error("Complex SVD could not complete the left singular basis");
                for (size_t row = 0; row < size; ++row)
                    column.left[row] = candidate[row] / norm;
            }
        }
    }

    void matrix<complex>::singularvalue(matrix<complex> &U, matrix<complex> &V,
            vector<double> &S) const
    {
        const size_t size = size_i();
        if (size_j() != size || U.size_i() != size || U.size_j() != size
                || V.size_i() != size || V.size_j() != size || S.size() != size || &U == &V)
            throw std::invalid_argument("Complex SVD requires square matrices of equal size and a matching vector");
        if (size == 0) return;

        std::vector<JacobiColumn> columns;
        columns.reserve(size);
        for (size_t index = 0; index < size; ++index) {
            columns.emplace_back(size);
            auto& column = columns.back();
            column.right[index] = 1.;
            for (size_t row = 0; row < size; ++row) {
                const complex entry = (*this)(row, index);
                if (!std::isfinite(entry.real()) || !std::isfinite(entry.imag()))
                    throw std::domain_error("Complex SVD requires finite entries");
                column.left[row] = WideComplex(entry.real(), entry.imag());
            }
            column.norm = columnNorm(column.left);
        }

        // Orthogonalize A V directly, retaining the relative accuracy of small columns.
        const long double tolerance = 8. * size * std::numeric_limits<long double>::epsilon();
        const size_t maxSweeps = std::max(size_t(30), 5 * size);
        bool converged = false;
        for (size_t sweep = 0; sweep < maxSweeps && !converged; ++sweep) {
            converged = true;
            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 1; j < size; ++j)
                    if (jacobiRotate(columns[i], columns[j], tolerance)) converged = false;
        }
        if (!converged) throw std::runtime_error("Complex SVD did not converge");

        for (auto& column : columns)
            if (column.norm != 0.)
                for (auto& entry : column.left) entry /= column.norm;
        completeLeftBasis(columns);

        std::vector<size_t> order(size);
        std::iota(order.begin(), order.end(), size_t(0));
        std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return columns[i].norm < columns[j].norm;
        });
        for (size_t index = 0; index < size; ++index) {
            const auto& column = columns[order[index]];
            S(index) = column.norm > std::numeric_limits<double>::max()
                    ? std::numeric_limits<double>::infinity() : static_cast<double>(column.norm);
            for (size_t row = 0; row < size; ++row) {
                U.assign(row, index, complex(column.left[row].real(), column.left[row].imag()));
                V.assign(row, index, complex(column.right[row].real(), column.right[row].imag()));
            }
        }
    }

    /** Conversion  */
    gsl_matrix_complex* matrix<complex>::as_gsl_type_ptr() const
    {
        return _matrix;
    }

    gsl_matrix_complex& matrix<complex>::as_gsl_type()
    {
        return const_cast<gsl_matrix_complex&> (*_matrix);
    }

    const gsl_matrix_complex& matrix<complex>::as_gsl_type() const
    {
        return const_cast<gsl_matrix_complex&> (*_matrix);
    }

    //  bool matrix<complex>::is_equal(const matrix<complex>& m1, const matrix<complex>& m2){
    //      return gsl_matrix_complex_equal(m1.as_gsl_type_ptr(),m2.as_gsl_type_ptr());
    //  }

    /** Unary minus */
    matrix<complex> matrix<complex>::operator-() const
    {
        matrix<complex> m1(_matrix);
        gsl_complex z1;
        GSL_SET_COMPLEX(&z1, -1., 0.);
        if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z1))
        {
            std::cout << "\n ** Error in matrix<complex> unary -" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Addition operator (matrix) */
    matrix<complex> matrix<complex>::operator+(const matrix<complex>& m) const
    {
        matrix<complex> m1(_matrix);
        if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> +" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Subtraction operator (matrix) */
    matrix<complex> matrix<complex>::operator-(const matrix<complex>& m) const
    {
        matrix<complex> m1(_matrix);
        matrix<complex> m2 = -m;
        if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> +" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Multiplication operator (matrix complex) */
    matrix<complex> matrix<complex>::operator*(const matrix<complex>& m) const
    {
        matrix<complex> m1(size_i(), m.size_j(), 0.);
        gsl_complex z1, z2;

        if (size_j() != m.size_i())
        {
            std::cout << "\n ** Size mismatch in matrix<complex> *" << std::endl;
            exit(EXIT_FAILURE);
        }
        GSL_SET_COMPLEX(&z1, 1., 0.);
        GSL_SET_COMPLEX(&z2, 0., 0.);
        if (gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, z1, _matrix,
                m.as_gsl_type_ptr(), z2, m1.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> *" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Addition operator (matrix) */
    matrix<complex> matrix<complex>::operator+(const matrix<double>& m) const
    {
        matrix<complex> m1(_matrix);
        matrix<complex> m2(m);
        if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> +" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Subtraction operator (matrix) */
    matrix<complex> matrix<complex>::operator-(const matrix<double>& m) const
    {
        matrix<complex> m1(_matrix);
        matrix<complex> m2(-m);
        if (gsl_matrix_complex_add(m1.as_gsl_type_ptr(), m2.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> +" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Multiplication operator (matrix double) */
    matrix<complex> matrix<complex>::operator*(const matrix<double>& m) const
    {
        matrix<complex> m1(m);
        return *this * m1;
    }

    /** Multiplication operator (vector complex) */
    vector<complex> matrix<complex>::operator*(const vector<complex>& v) const
    {
        gsl_complex z1, z2;
        vector<complex> v1(size_i(), 0.);

        if (size_j() != v.size())
        {
            std::cout << "\n ** Size mismatch in matrix<complex> * (vector complex)"
                    << std::endl;
            exit(EXIT_FAILURE);
        }
        GSL_SET_COMPLEX(&z1, 1., 0.);
        GSL_SET_COMPLEX(&z2, 0., 0.);
        if (gsl_blas_zgemv(CblasNoTrans, z1, _matrix, v.as_gsl_type_ptr(),
                z2, v1.as_gsl_type_ptr()))
        {
            std::cout << "\n ** Error in matrix<complex> * (vector complex)"
                    << std::endl;
            exit(EXIT_FAILURE);
        }
        return v1;
    }

    /** Multiplication operator (vector double) */
    vector<complex> matrix<complex>::operator*(const vector<double>& v) const
    {
        vector<complex> v1(v);
        return *this * v1;
    }

    /** Addition assignment (matrix) */
    matrix<complex>& matrix<complex>::operator+=(const matrix<complex>& m)
    {
        *this = *this +m;
        return *this;
    }

    /** Subtraction assignment (matrix) */
    matrix<complex>& matrix<complex>::operator-=(const matrix<complex>& m)
    {
        *this = *this -m;
        return *this;
    }

    /** Multiplication assignment (matrix) */
    matrix<complex>& matrix<complex>::operator*=(const matrix<complex>& m)
    {
        if (!((size_i() == size_j()) && (size_i() == m.size_i()) &&
                (size_j() == m.size_j())))
        {
            std::cout << "\n ** Size mismatch in matrix<complex> *= (matrix)"
                    << std::endl;
            exit(EXIT_FAILURE);
        }
        *this = *this * m;
        return *this;
    }

    /** Addition operator (complex) */
    matrix<complex> matrix<complex>::operator+(const complex& z) const
    {
        matrix<complex> m1(_matrix);
        if (gsl_matrix_complex_add_constant(m1.as_gsl_type_ptr(), z.as_gsl_type()))
        {
            std::cout << "\n ** Error in matrix<complex> + (complex)" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Subtraction operator (complex) */
    matrix<complex> matrix<complex>::operator-(const complex& z) const
    {
        matrix<complex> m1(_matrix);
        if (gsl_matrix_complex_add_constant(m1.as_gsl_type_ptr(), (-z).as_gsl_type()))
        {
            std::cout << "\n ** Error in matrix<complex> - (complex)" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Multiplication operator (complex) */
    matrix<complex> matrix<complex>::operator*(const complex& z) const
    {
        matrix<complex> m1(_matrix);
        if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z.as_gsl_type()))
        {
            std::cout << "\n ** Error in matrix<complex> * (complex)" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Division operator (complex) */
    matrix<complex> matrix<complex>::operator/(const complex& z) const
    {
        matrix<complex> m1(_matrix);
        if (gsl_matrix_complex_scale(m1.as_gsl_type_ptr(), z.inverse().as_gsl_type()))
        {
            std::cout << "\n ** Error in matrix<complex> / (complex)" << std::endl;
            exit(EXIT_FAILURE);
        }
        return m1;
    }

    /** Addition assignment (complex) */
    matrix<complex>& matrix<complex>::operator+=(const complex& z)
    {
        *this = *this +z;
        return *this;
    }

    /** Subtraction assignment (complex) */
    matrix<complex>& matrix<complex>::operator-=(const complex& z)
    {
        *this = *this -z;
        return *this;
    }

    /** Multiplication assignment (complex) */
    matrix<complex>& matrix<complex>::operator*=(const complex& z)
    {
        *this = *this * z;
        return *this;
    }

    /** Division assignment (complex) */
    matrix<complex>& matrix<complex>::operator/=(const complex& z)
    {
        *this = *this / z;
        return *this;
    }

    /** Addition operator (double) */
    matrix<complex> matrix<complex>::operator+(const double& a) const
    {
        complex z(a);
        return *this +z;
    }

    /** Subtraction assignment (double) */
    matrix<complex> matrix<complex>::operator-(const double& a) const
    {
        complex z(a);
        return *this -z;
    }

    /** Multiplication operator (double) */
    matrix<complex> matrix<complex>::operator*(const double& a) const
    {
        complex z(a);
        return *this * z;
    }

    /** Division operator (double) */
    matrix<complex> matrix<complex>::operator/(const double& a) const
    {
        complex z(a);
        return *this / z;
    }

    /** Addition assignment (double) */
    matrix<complex>& matrix<complex>::operator+=(const double& a)
    {
        *this = *this +a;
        return *this;
    }

    /** Subtraction assignment (double) */
    matrix<complex>& matrix<complex>::operator-=(const double& a)
    {
        *this = *this -a;
        return *this;
    }

    /** Multiplication assignment (double) */
    matrix<complex>& matrix<complex>::operator*=(const double& a)
    {
        *this = *this * a;
        return *this;
    }

    /** Division assignment (double) */
    matrix<complex>& matrix<complex>::operator/=(const double& a)
    {
        *this = *this / a;
        return *this;
    }

    /** Comparison == (matrix) */
    bool matrix<complex>::operator==(const matrix<complex>& a) const
    {
        if (a.size_i() != size_i() || a.size_j() != size_j())
        {
            std::cout << "\n Error in matrix<complex>::operator== (matrix): cannot compare matrices of different size" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i < size_i(); i++)
            for (size_t j = 0; j < size_j(); j++)
                if (a(i, j) != (*this)(i, j)) return (false);
        return (true);
    }

    /** friend functions */
    std::ostream& operator<<(std::ostream& output, const matrix<complex>& m)
    {
        size_t i, j;
        for (i = 0; i < m.size_i(); i++)
        {
            output << std::endl;
            output << "\t(";
            for (j = 0; j < m.size_j() - 1; j++)
                output << m(i, j) << ",";
            output << m(i, j) << ")";
        }
        return output;
    }

    matrix<complex> operator+(const matrix<double> m1, const matrix<complex> m2)
    {
        return m2 + m1;
    }

    matrix<complex> operator-(const matrix<double> m1, const matrix<complex> m2)
    {
        return -m2 + m1;
    }

    matrix<complex> operator+(const complex& z, const matrix<complex> m)
    {
        return m + z;
    }

    matrix<complex> operator-(const complex& z, const matrix<complex> m)
    {
        return -m + z;
    }

    matrix<complex> operator*(const complex& z, const matrix<complex> m)
    {
        return m * z;
    }

    vector<complex> operator*(const vector<complex>& v, const matrix<complex> m)
    {
        return m.transpose() * v;
    }

    vector<complex> operator*(const vector<double>& v, const matrix<complex> m)
    {
        return m.transpose() * v;
    }

    matrix<complex> operator+(const double& a, const matrix<complex> m)
    {
        return m + a;
    }

    matrix<complex> operator-(const double& a, const matrix<complex> m)
    {
        return -m + a;
    }

    matrix<complex> operator*(const double& a, const matrix<complex> m)
    {
        return m * a;
    }
}
