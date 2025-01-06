#pragma once

#include "SO2.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <math.h>
#include <limits>

using namespace Eigen;

/**
 * @brief Class representing a member of the \f$SE(2)\f$ manifold, or a 2D rigid body transform.
 */
template<typename T>
class SE2
{
private:
    typedef Matrix<T, 1, 1> Vec1T;
    typedef Matrix<T, 2, 1> Vec2T;
    typedef Matrix<T, 3, 1> Vec3T;
    typedef Matrix<T, 4, 1> Vec4T;
    typedef Matrix<T, 2, 2> Mat2T;
    typedef Matrix<T, 3, 3> Mat3T;

    T buf_[4];

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * @brief Memory-mapped array representing all transform fields in \f$\begin{bmatrix}\boldsymbol{t} &
     * \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    Map<Vec4T> arr_;
    /**
     * @brief Memory-mapped array representing only the translation component of the transform \f$\boldsymbol{t}\f$.
     */
    Map<Vec2T> t_;
    /**
     * @brief Memory-mapped array representing only the rotation component of the transform \f$\boldsymbol{q}\f$.
     */
    SO2<T> q_;

    /**
     * @brief Obtain a random \f$SE(2)\f$ transform.
     *
     * \f$\boldsymbol{q}\f$ is normalized and \f$0\leq t_i \leq 1\f$.
     */
    static SE2 random()
    {
        SE2 x;
        x.arr_.setRandom();
        x.q_.normalize();
        return x;
    }

    /**
     * @brief Obtain an identity \f$SE(2)\f$ transform.
     */
    static SE2 identity()
    {
        SE2 x;
        x.arr_.setZero();
        x.arr_(2) = (T)1.0;
        return x;
    }

    /**
     * @brief Obtain a transform full of NaNs.
     */
    static SE2 nans()
    {
        SE2 x;
        x.arr_.setConstant(std::numeric_limits<T>::quiet_NaN());
        return x;
    }

    /**
     * @brief Convert a transform matrix \f$\begin{bmatrix}\boldsymbol{R} & \boldsymbol{t} \\ \boldsymbol{0} &
     * 1\end{bmatrix}\f$ into a \f$SE(2)\f$ transform.
     */
    static SE2 fromH(const Mat3T& m)
    {
        return SE2::fromVecAndRot(m.template block<2, 1>(0, 2), SO2<T>::fromR(m.template block<2, 2>(0, 0)));
    }

    /**
     * @brief Construct a transform from the individual fields.
     * @param tx First component of the translation \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     * @param ty Second component of the translation \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     * @param qw Real component of the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     * @param qx Imaginary component of the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     */
    static SE2 fromVecAndRot(const T tx, const T ty, const T qw, const T qx)
    {
        SE2 x;
        x.arr_ << tx, ty, qw, qx;
        return x;
    }

    /**
     * @brief Construct a transform from a translation vector and rotation fields vector.
     * @param tvec The translation vector \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     * @param qvec The rotation represented as an array \f$\boldsymbol{q}\in\mathbb{R}^2\f$.
     */
    static SE2 fromVecAndRot(const Vec2T& tvec, const Vec2T& qvec)
    {
        SE2 x;
        x.arr_ << tvec(0), tvec(1), qvec(0), qvec(1);
        return x;
    }

    /**
     * @brief Construct a transform from a translation vector and a rotation object.
     * @param t The translation vector \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     * @param q The rotation \f$\boldsymbol{q}\in SO(2)\f$.
     */
    static SE2 fromVecAndRot(const Vec2T& t, const SO2<T>& q)
    {
        SE2 x;
        x.arr_ << t(0), t(1), q.w(), q.x();
        return x;
    }

    /**
     * @brief Create a transform (with garbage data).
     */
    SE2() : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 2) {}

    /**
     * @brief Create a transform from an array representing all transform fields in \f$\begin{bmatrix}\boldsymbol{t} &
     * \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    SE2(const Ref<const Vec4T>& arr) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 2)
    {
        arr_ = arr;
    }

    /**
     * @brief Copy constructor from another transform.
     */
    SE2(const SE2& x) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 2)
    {
        arr_ = x.arr_;
    }

    /**
     * @brief Create a transform from a pointer array representing all transform fields in
     * \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    SE2(const T* data) : arr_(const_cast<T*>(data)), t_(arr_.data()), q_(arr_.data() + 2) {}

    /**
     * @brief Access a field from \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline T& operator[](int i)
    {
        return arr_[i];
    }

    /**
     * @brief Access the translation vector \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     */
    inline const Map<Vec2T>& t() const
    {
        return t_;
    }

    /**
     * @brief Access the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     */
    inline const SO2<T>& q() const
    {
        return q_;
    }

    /**
     * @brief Access the translation vector \f$\boldsymbol{t}\in\mathbb{R}^2\f$.
     */
    inline Map<Vec2T>& t()
    {
        return t_;
    }

    /**
     * @brief Access the rotation \f$\boldsymbol{q}\in SO(2)\f$.
     */
    inline SO2<T>& q()
    {
        return q_;
    }

    /**
     * @brief Access all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline const Vec4T elements() const
    {
        return arr_;
    }

    /**
     * @brief Access all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline Vec4T array() const
    {
        return arr_;
    }

    /**
     * @brief Access pointer to all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline T* data()
    {
        return arr_.data();
    }

    /**
     * @brief Access pointer to all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline const T* data() const
    {
        return arr_.data();
    }

    /**
     * @brief Get a deep copy of the current transform.
     */
    SE2 copy() const
    {
        SE2 tmp;
        tmp.arr_ = arr_;
        return tmp;
    }

    /**
     * @brief Convert the transform to matrix representation \f$\begin{bmatrix}\boldsymbol{R} & \boldsymbol{t}
     * \\ \boldsymbol{0} & 1\end{bmatrix}\in\mathbb{R}^{3\times 3}\f$
     */
    Mat3T H() const
    {
        Mat3T out;
        out << q_.R(), t_, (T)0., (T)0., (T)1.;
        return out;
    }

    /**
     * @brief Obtain the inverse transform \f$\boldsymbol{T}_A^B\rightarrow \boldsymbol{T}_B^A\f$.
     */
    SE2 inverse() const
    {
        SE2    x;
        SO2<T> q_inv = q_.inverse();
        return SE2::fromVecAndRot(-(q_inv * t_), q_inv);
    }

    /**
     * @brief Invert the current transform \f$\boldsymbol{T}_A^B\rightarrow \boldsymbol{T}_B^A\f$.
     */
    SE2& invert()
    {
        q().invert();
        t() = -(q() * t());
        return *this;
    }

    /**
     * @brief Implementation of group composition: \f$\boldsymbol{T}_B^C \otimes \boldsymbol{T}_A^B\rightarrow
     * \boldsymbol{T}_A^C\f$.
     */
    template<typename Tout = T, typename T2>
    SE2<Tout> otimes(const SE2<T2>& x) const
    {
        SE2<Tout> xout;
        xout.t() = t_ + q_ * x.t_;
        xout.q() = q_ * x.q_;
        return xout;
    }

    /**
     * @brief Implementation of tangent space group perturbations: \f$\boldsymbol{T}_A^B\oplus \boldsymbol{t}_B^{B'}
     * \rightarrow \boldsymbol{T}_A^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    SE2<Tout> oplus(const Matrix<T2, 3, 1>& delta) const
    {
        return otimes<Tout, T2>(SE2<T2>::Exp(delta));
    }

    /**
     * @brief Implementation of group subtraction: \f$\boldsymbol{T}_A^B\ominus \boldsymbol{T}_A^{B'} \rightarrow
     * \boldsymbol{t}_B^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 3, 1> ominus(const SE2<T2>& x) const
    {
        return SE2<Tout>::Log(x.inverse().template otimes<Tout>(*this));
    }

    /**
     * @brief Copy constructor.
     */
    SE2& operator=(const SE2& x)
    {
        arr_ = x.elements();
        return *this;
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SE2 operator*(const SE2<T2>& x) const
    {
        return otimes(x);
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SE2& operator*=(const SE2<T2>& x)
    {
        arr_ = otimes(x).elements();
        return *this;
    }

    /**
     * @brief Scale a transform by a scalar.
     *
     * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a transform.
     */
    SE2& operator*=(const double& s)
    {
        arr_ = SE2::Exp(s * SE2::Log(*this)).elements();
        return *this;
    }

    /**
     * @brief Scale a transform by a scalar.
     *
     * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a transform.
     */
    SE2& operator/=(const double& s)
    {
        arr_ = SE2::Exp(SE2::Log(*this) / s).elements();
        return *this;
    }

    /**
     * @brief Scale a transform by a scalar.
     *
     * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a transform.
     */
    SE2 operator/(const double& s) const
    {
        SE2 qs;
        qs.arr_ = SE2::Exp(SE2::Log(*this) / s).elements();
        return qs;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{T}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 2, 1> operator*(const Matrix<T2, 2, 1>& v) const
    {
        return q_ * v + t_;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{T}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    Vec2T operator*(const Vec2T& v) const
    {
        return q_ * v + t_;
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SE2 operator+(const Vec3T& v) const
    {
        return oplus(v);
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SE2& operator+=(const Vec3T& v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }

    /**
     * @brief Invocation of ominus via subtraction.
     */
    template<typename T2>
    Vec3T operator-(const SE2<T2>& x) const
    {
        return ominus(x);
    }

    /**
     * @brief Hat operator implementation, which coverts the tangent-space vector representation to the corresponding
     * Lie algebra: \f$\mathbb{R}^3\rightarrow \mathfrak{se}(2)\f$.
     */
    static Mat3T hat(const Vec3T& omega)
    {
        Mat3T Omega;
        Omega << SO2<T>::hat(omega.template block<1, 1>(2, 0)), omega.template block<2, 1>(0, 0), (T)0., (T)0., (T)0.;
        return Omega;
    }

    /**
     * @brief Vee operator implementation, which coverts the Lie algebra representation to a tangent-space vector
     * representation: \f$\mathfrak{se}(2) \rightarrow \mathbb{R}^3\f$.
     */
    static Vec3T vee(const Mat3T& Omega)
    {
        Vec3T omega;
        omega << Omega.template block<2, 1>(0, 2), SO2<T>::vee(Omega.template block<2, 2>(0, 0));
        return omega;
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SE(2) \rightarrow \mathfrak{se}(2)\f$.
     */
    static Mat3T log(const SE2& x)
    {
        return SE2::hat(SE2::Log(x));
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SE(2) \rightarrow \mathbb{R}^3\f$.
     */
    static Vec3T Log(const SE2& x)
    {
        Vec1T w  = SO2<T>::Log(x.q_);
        T     th = w.norm();

        Mat2T leftJacobianInverse;
        if (th > (T)1e-4)
        {
            T     a = sin(th) / ((T)1. - cos(th));
            Mat2T skew1;
            skew1 << (T)0, (T)-1., (T)1., (T)0;
            leftJacobianInverse = th / (T)2. * (a * Mat2T::Identity() - skew1);
        }
        else
        {
            leftJacobianInverse = Mat2T::Identity();
        }

        Vec3T out;
        out << leftJacobianInverse * x.t_, w;

        return out;
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathfrak{se}(2) \rightarrow SE(2)\f$.
     */
    static SE2 exp(const Mat3T& Omega)
    {
        return SE2::Exp(SE2::vee(Omega));
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathbb{R}^3 \rightarrow SE(2)\f$.
     */
    static SE2 Exp(const Vec3T& omega)
    {
        Vec2T  rho = omega.template block<2, 1>(0, 0);
        Vec1T  w   = omega.template block<1, 1>(2, 0);
        SO2<T> q   = SO2<T>::Exp(w);
        T      th  = w.norm();

        Mat2T leftJacobian;
        if (abs(th) > (T)1e-4)
        {
            T     a = sin(th) / th;
            T     b = ((T)1. - cos(th)) / th;
            Mat2T skew1;
            skew1 << (T)0, (T)-1., (T)1., (T)0;
            leftJacobian = a * Mat2T::Identity() + b * skew1;
        }
        else
        {
            leftJacobian = Mat2T::Identity();
        }

        return SE2::fromVecAndRot(leftJacobian * rho, q);
    }

    /**
     * @brief Cast the underlying numeric type.
     */
    template<typename T2>
    SE2<T2> cast() const
    {
        SE2<T2> x;
        x.arr_ = arr_.template cast<T2>();
        return x;
    }
};

/**
 * @brief Scale a transform by a scalar.
 *
 * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a transform.
 */
template<typename T>
SE2<T> operator*(const double& l, const SE2<T>& r)
{
    SE2<T> lr;
    lr.arr_ = SE2<T>::Exp(l * SE2<T>::Log(r)).elements();
    return lr;
}

/**
 * @brief Scale a transform by a scalar.
 *
 * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a transform.
 */
template<typename T>
SE2<T> operator*(const SE2<T>& l, const double& r)
{
    SE2<T> lr;
    lr.arr_ = SE2<T>::Exp(r * SE2<T>::Log(l)).elements();
    return lr;
}

/**
 * @brief Render the transform in a stream.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SE2<T>& x)
{
    os << "SE(2): [ " << x.t_.x() << "i, " << x.t_.y() << "j ] [ " << x.q_.w() << ", " << x.q_.x() << "i ]";
    return os;
}

typedef SE2<double> SE2d;
