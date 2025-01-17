#pragma once

#include <Eigen/Core>
#include <iostream>
#include <math.h>

using namespace Eigen;

/**
 * @brief Class representing a member of the \f$SO(2)\f$ manifold, or a 2D rotation.
 */
template<typename T>
class SO2
{
private:
    typedef Matrix<T, 1, 1> Vec1T;
    typedef Matrix<T, 2, 1> Vec2T;
    typedef Matrix<T, 2, 2> Mat2T;

    T buf_[2];

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * @brief Memory-mapped array representing all rotation fields in \f$\boldsymbol{q}\f$.
     */
    Map<Vec2T> arr_;

    /**
     * @brief Obtain a random rotation.
     * @return A random rotation \f$\boldsymbol{q}_B^W\in SO(2)\f$.
     *
     *  The rotation \f$\mathbf{q}_B^W\f$ will be normalized.
     */
    static SO2 random()
    {
        SO2 q;
        q.arr_.setRandom();
        q.normalize();
        return q;
    }

    /**
     * @brief Obtain an identity \f$SO(2)\f$ rotation.
     */
    static SO2 identity()
    {
        SO2 q;
        q.arr_ << 1., 0.;
        return q;
    }

    /**
     * @brief Obtain a rotation full of NaNs.
     */
    static SO2 nans()
    {
        SO2 x;
        x.arr_.setConstant(std::numeric_limits<T>::quiet_NaN());
        return x;
    }

    /**
     * @brief Convert an angle (in radians) into a rotation.
     */
    static SO2 fromAngle(const T angle)
    {
        SO2 q;
        q.arr_ << cos(angle), sin(angle);
        return q;
    }

    /**
     * @brief Convert a rotation matrix \f$\boldsymbol{R}\in\mathbb{R}^{2\times 2}\f$ to a \f$SO(2)\f$ rotation.
     */
    static SO2 fromR(const Mat2T& m)
    {
        SO2 q;
        q.arr_ << m(0, 0), m(1, 0);
        return q;
    }

    /**
     * @brief Given two unit vectors \f$\boldsymbol{u},\boldsymbol{v}\in \mathbb{R}^2\f$, returns the rotation that
     * rotates \f$\boldsymbol{u}\rightarrow\boldsymbol{v}\f$.
     */
    static SO2 fromTwoUnitVectors(const Vec2T& u, const Vec2T& v)
    {
        SO2 q;
        T   d = u.dot(v);

        if (d < (T)0.99999999 && d > (T)-0.99999999)
        {
            q.arr_ << d, u(0) * v(1) - u(1) * v(0);
        }
        else if (d < (T)-0.99999999)
        {
            q.arr_ << -1., 0.;
        }
        else
        {
            q = SO2::identity();
        }

        return q;
    }

    /**
     * @brief Construct a rotation from the individual fields.
     * @param qw Real component of the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     * @param qx Imaginary component of the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     */
    static SO2 fromComplex(const T qw, const T qx)
    {
        SO2 q;
        q.arr_ << qw, qx;
        return q;
    }

    /**
     * @brief Construct a rotation from a rotation fields vector.
     * @param qvec The rotation represented as an array \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    static SO2 fromComplex(const Vec2T& qvec)
    {
        SO2 q;
        q.arr_ << qvec(0), qvec(1);
        return q;
    }

    /**
     * @brief Create a rotation (with garbage data).
     */
    SO2() : arr_(buf_) {}

    /**
     * @brief Create a transform from an array representing all rotation fields in \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    SO2(const Ref<const Vec2T>& arr) : arr_(buf_)
    {
        arr_ = arr;
    }

    /**
     * @brief Copy constructor from another rotation.
     */
    SO2(const SO2& q) : arr_(buf_)
    {
        arr_ = q.arr_;
    }

    /**
     * @brief Create a rotation from a pointer array representing all rotation fields in
     * \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    SO2(const T* data) : arr_(const_cast<T*>(data)) {}

    /**
     * @brief Access a field from \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    inline T& operator[](int i)
    {
        return arr_[i];
    }

    /**
     * @brief Access the real element of the rotation.
     */
    inline const T& w() const
    {
        return arr_(0);
    }

    /**
     * @brief Access the imaginary element of the rotation.
     */
    inline const T& x() const
    {
        return arr_(1);
    }

    /**
     * @brief Access the real element of the rotation.
     */
    inline T& w()
    {
        return arr_(0);
    }

    /**
     * @brief Access the imaginary element of the rotation.
     */
    inline T& x()
    {
        return arr_(1);
    }

    /**
     * @brief Access all elements of \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    inline const Vec2T elements() const
    {
        return arr_;
    }

    /**
     * @brief Access all elements of \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    inline Vec2T array() const
    {
        return arr_;
    }

    /**
     * @brief Access pointer to all elements of \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    inline T* data()
    {
        return arr_.data();
    }

    /**
     * @brief Access pointer to all elements of \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    inline const T* data() const
    {
        return arr_.data();
    }

    /**
     * @brief Get a deep copy of the current rotation.
     */
    SO2 copy() const
    {
        SO2 tmp;
        tmp.arr_ = arr_;
        return tmp;
    }

    /**
     * @brief Normalize the elements of \f$\boldsymbol{q}\f$ to make it a valid rotation.
     */
    void normalize()
    {
        arr_ /= arr_.norm();
    }

    /**
     * @brief Obtain a normalized copy of the current rotation.
     */
    SO2 normalized()
    {
        SO2 tmp = copy();
        tmp.normalize();
        return tmp;
    }

    /**
     * @brief Convert the rotation to matrix representation \f$\begin{bmatrix}\cos\theta & -\sin\theta
     * \\ \sin\theta & \cos\theta\end{bmatrix}\f$.
     */
    Mat2T R() const
    {
        T     costh = w();
        T     sinth = x();
        Mat2T out;
        out << costh, -sinth, sinth, costh;
        return out;
    }

    /**
     * @brief Obtain the inverse rotation \f$\boldsymbol{R}_A^B\rightarrow \boldsymbol{R}_B^A\f$.
     */
    SO2 inverse() const
    {
        SO2 q;
        q.arr_(0) = arr_(0);
        q.arr_(1) = -arr_(1);
        return q;
    }

    /**
     * @brief Invert the current rotation \f$\boldsymbol{R}_A^B\rightarrow \boldsymbol{R}_B^A\f$.
     */
    SO2& invert()
    {
        arr_(1) *= (T)-1.0;
        return *this;
    }

    /**
     * @brief Obtain the equivalent scalar angle \f$\theta\f$ of the rotation.
     *
     * \f$-\pi \leq \theta \leq \pi\f$.
     */
    T angle() const
    {
        return atan2(x(), w());
    }

    /**
     * @brief Implementation of group composition: \f$\boldsymbol{q}_B^C \otimes \boldsymbol{q}_A^B\rightarrow
     * \boldsymbol{q}_A^C\f$.
     */
    template<typename Tout = T, typename T2>
    SO2<Tout> otimes(const SO2<T2>& q) const
    {
        SO2<Tout> qout;
        qout.arr_ << w() * q.w() - x() * q.x(), w() * q.x() + x() * q.w();
        return qout;
    }

    /**
     * @brief Implementation of tangent space group perturbations: \f$\boldsymbol{q}_A^B\oplus
     * \boldsymbol{\theta}_B^{B'} \rightarrow \boldsymbol{q}_A^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    SO2<Tout> oplus(const Matrix<T2, 1, 1>& delta) const
    {
        return otimes<Tout, T2>(SO2<T2>::Exp(delta));
    }

    /**
     * @brief Implementation of group subtraction: \f$\boldsymbol{q}_A^B\ominus \boldsymbol{q}_A^{B'} \rightarrow
     * \boldsymbol{\theta}_B^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 1, 1> ominus(const SO2<T2>& q) const
    {
        SO2<Tout> dq = q.inverse().template otimes<Tout>(*this);
        return SO2<Tout>::Log(dq);
    }

    /**
     * @brief Copy constructor.
     */
    SO2& operator=(const SO2& q)
    {
        arr_ = q.elements();
        return *this;
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SO2 operator*(const SO2<T2>& q) const
    {
        return otimes(q);
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SO2& operator*=(const SO2<T2>& q)
    {
        arr_ = otimes(q).elements();
        return *this;
    }

    /**
     * @brief Scale a rotation by a scalar.
     *
     * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a rotation.
     */
    SO2& operator*=(const double& s)
    {
        arr_ = SO2::Exp(s * SO2::Log(*this)).elements();
        return *this;
    }

    /**
     * @brief Scale a rotation by a scalar.
     *
     * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a rotation.
     */
    SO2& operator/=(const double& s)
    {
        arr_ = SO2::Exp(SO2::Log(*this) / s).elements();
        return *this;
    }

    /**
     * @brief Scale a rotation by a scalar.
     *
     * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a rotation.
     */
    SO2 operator/(const double& s) const
    {
        SO2 qs;
        qs.arr_ = SO2::Exp(SO2::Log(*this) / s).elements();
        return qs;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{q}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 2, 1> operator*(const Matrix<T2, 2, 1>& v) const
    {
        Vec2T out;
        out << w() * v.x() - x() * v.y(), w() * v.y() + x() * v.x();
        return out;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{q}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    Vec2T operator*(const Vec2T& v) const
    {
        Vec2T out;
        out << w() * v.x() - x() * v.y(), w() * v.y() + x() * v.x();
        return out;
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SO2 operator+(const Vec1T& v) const
    {
        return oplus(v);
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SO2& operator+=(const Vec1T& v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }

    /**
     * @brief Invocation of ominus via subtraction.
     */
    template<typename T2>
    Vec1T operator-(const SO2<T2>& q) const
    {
        return ominus(q);
    }

    /**
     * @brief Hat operator implementation, which coverts the tangent-space vector representation to the corresponding
     * Lie algebra: \f$\mathbb{R}\rightarrow \mathfrak{so}(2)\f$.
     */
    static Mat2T hat(const Vec1T& omega)
    {
        Mat2T Omega;
        Omega << (T)0., -omega.x(), omega.x(), (T)0.;
        return Omega;
    }

    /**
     * @brief Vee operator implementation, which coverts the Lie algebra representation to a tangent-space vector
     * representation: \f$\mathfrak{so}(2) \rightarrow \mathbb{R}\f$.
     */
    static Vec1T vee(const Mat2T& Omega)
    {
        Vec1T omega;
        omega << Omega(1, 0);
        return omega;
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SO(2) \rightarrow \mathfrak{so}(2)\f$.
     */
    static Mat2T log(const SO2& q)
    {
        return hat(SO2::Log(q));
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SO(2) \rightarrow \mathbb{R}\f$.
     */
    static Vec1T Log(const SO2& q)
    {
        Vec1T out;
        out << q.angle();
        return out;
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathfrak{so}(2) \rightarrow SO(2)\f$.
     */
    static SO2 exp(const Mat2T& Omega)
    {
        return SO2::Exp(vee(Omega));
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathbb{R} \rightarrow SO(2)\f$.
     */
    static SO2 Exp(const Vec1T& omega)
    {
        return SO2::fromAngle(omega.x());
    }

    /**
     * @brief Cast the underlying numeric type.
     */
    template<typename T2>
    SO2<T2> cast() const
    {
        SO2<T2> q;
        q.arr_ = arr_.template cast<T2>();
        return q;
    }
};

/**
 * @brief Scale a rotation by a scalar.
 *
 * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a rotation.
 */
template<typename T>
SO2<T> operator*(const double& l, const SO2<T>& r)
{
    SO2<T> lr;
    lr.arr_ = SO2<T>::Exp(l * SO2<T>::Log(r)).elements();
    return lr;
}

/**
 * @brief Scale a rotation by a scalar.
 *
 * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a rotation.
 */
template<typename T>
SO2<T> operator*(const SO2<T>& l, const double& r)
{
    SO2<T> lr;
    lr.arr_ = SO2<T>::Exp(r * SO2<T>::Log(l)).elements();
    return lr;
}

/**
 * @brief Render the rotation in a stream.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SO2<T>& q)
{
    os << "SO(2): [ " << q.w() << ", " << q.x() << "i ]";
    return os;
}

typedef SO2<double> SO2d;
