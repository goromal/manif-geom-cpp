#pragma once

#include "SO3.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <math.h>

using namespace Eigen;

/**
 * @brief Class representing a member of the \f$SE(3)\f$ manifold, or a 3D rigid body transform.
 */
template<typename T>
class SE3
{
private:
    typedef Matrix<T, 3, 1> Vec3T;
    typedef Matrix<T, 4, 1> Vec4T;
    typedef Matrix<T, 6, 1> Vec6T;
    typedef Matrix<T, 7, 1> Vec7T;
    typedef Matrix<T, 3, 3> Mat3T;
    typedef Matrix<T, 4, 4> Mat4T;
    typedef Matrix<T, 6, 6> Mat6T;

    T buf_[7];

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * @brief Memory-mapped array representing all transform fields in \f$\begin{bmatrix}\boldsymbol{t} &
     * \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    Map<Vec7T> arr_;
    /**
     * @brief Memory-mapped array representing only the translation component of the transform \f$\boldsymbol{t}\f$.
     */
    Map<Vec3T> t_;
    /**
     * @brief Memory-mapped array representing only the rotation component of the transform \f$\boldsymbol{q}\f$.
     */
    SO3<T> q_;

    /**
     * @brief Obtain a random rigid body transform.
     * @return A random rigid body transform \f$\mathbf{X}_B^W\in SE(3)\f$.
     *
     *  The rotation component \f$\mathbf{q}_B^W\f$ will be normalized, but the translation component
     * \f$\mathbf{t}_{B/W}^W\f$ will not (each component will be between 0 and 1).
     */
    static SE3 random()
    {
        SE3 x;
        x.arr_.setRandom();
        x.q_.normalize();
        return x;
    }

    /**
     * @brief Obtain an identity \f$SE(3)\f$ transform.
     */
    static SE3 identity()
    {
        SE3 x;
        x.arr_.setZero();
        x.arr_(3) = (T)1.0;
        return x;
    }

    /**
     * @brief Obtain a transform full of NaNs.
     */
    static SE3 nans()
    {
        SE3 x;
        x.arr_.setConstant(std::numeric_limits<T>::quiet_NaN());
        return x;
    }

    /**
     * @brief Convert a transform matrix \f$\begin{bmatrix}\boldsymbol{R} & \boldsymbol{t} \\ \boldsymbol{0} &
     * 1\end{bmatrix}\f$ into a \f$SE(3)\f$ transform.
     */
    static SE3 fromH(const Mat4T& m)
    {
        return SE3::fromVecAndQuat(m.template block<3, 1>(0, 3), SO3<T>::fromR(m.template block<3, 3>(0, 0)));
    }

    /**
     * @brief Construct a transform from the individual fields.
     * @param tx First component of the translation \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param ty Second component of the translation \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param tz Third component of the translation \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param qw Real component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qx First imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qy Second imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qz Third imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     */
    static SE3 fromVecAndQuat(const T tx, const T ty, const T tz, const T qw, const T qx, const T qy, const T qz)
    {
        SE3 x;
        x.arr_ << tx, ty, tz, qw, qx, qy, qz;
        return x;
    }

    /**
     * @brief Construct a transform from a translation vector and rotation fields vector.
     * @param tvec The translation vector \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param qvec The rotation represented as an array \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    static SE3 fromVecAndQuat(const Vec3T& tvec, const Vec4T& qvec)
    {
        SE3 x;
        x.arr_ << tvec(0), tvec(1), tvec(2), qvec(0), qvec(1), qvec(2), qvec(3);
        return x;
    }

    /**
     * @brief Construct a transform from a translation vector and a rotation object.
     * @param t The translation vector \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param q The rotation \f$\boldsymbol{q}\in SO(3)\f$.
     */
    static SE3 fromVecAndQuat(const Vec3T& t, const SO3<T>& q)
    {
        SE3 x;
        x.arr_ << t(0), t(1), t(2), q.w(), q.x(), q.y(), q.z();
        return x;
    }

    /**
     * @brief Construct a transform from a translation vector and an Eigen Quaternion object.
     * @param t The translation vector \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     * @param q The Eigen Quaternion \f$\boldsymbol{q}\in SO(3)\f$.
     */
    static SE3 fromVecAndQuat(const Vec3T& tvec, const Quaternion<T> quat)
    {
        SE3 x;
        x.arr_ << tvec(0), tvec(1), tvec(2), quat.w(), quat.x(), quat.y(), quat.z();
        return x;
    }

    /**
     * @brief Create a transform (with garbage data).
     */
    SE3() : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3) {}

    /**
     * @brief Create a transform from an array representing all transform fields in \f$\begin{bmatrix}\boldsymbol{t} &
     * \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    SE3(const Ref<const Vec7T>& arr) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3)
    {
        arr_ = arr;
    }

    /**
     * @brief Copy constructor from another transform.
     */
    SE3(const SE3& x) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3)
    {
        arr_ = x.arr_;
    }

    /**
     * @brief Create a transform from a pointer array representing all transform fields in
     * \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    SE3(const T* data) : arr_(const_cast<T*>(data)), t_(arr_.data()), q_(arr_.data() + 3) {}

    /**
     * @brief Access a field from \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline T& operator[](int i)
    {
        return arr_[i];
    }

    /**
     * @brief Access the translation vector \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     */
    inline const Map<Vec3T>& t() const
    {
        return t_;
    }

    /**
     * @brief Access the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     */
    inline const SO3<T>& q() const
    {
        return q_;
    }

    /**
     * @brief Access the translation vector \f$\boldsymbol{t}\in\mathbb{R}^3\f$.
     */
    inline Map<Vec3T>& t()
    {
        return t_;
    }

    /**
     * @brief Access the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     */
    inline SO3<T>& q()
    {
        return q_;
    }

    /**
     * @brief Access all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline const Vec7T elements() const
    {
        return arr_;
    }

    /**
     * @brief Access all elements of \f$\begin{bmatrix}\boldsymbol{t} & \boldsymbol{q}\end{bmatrix}^\top\f$.
     */
    inline Vec7T array() const
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
    SE3 copy() const
    {
        SE3 tmp;
        tmp.arr_ = arr_;
        return tmp;
    }

    /**
     * @brief Convert the transform to matrix representation \f$\begin{bmatrix}\boldsymbol{R} & \boldsymbol{t}
     * \\ \boldsymbol{0} & 1\end{bmatrix}\in\mathbb{R}^{4\times 4}\f$
     */
    Mat4T H() const
    {
        Mat4T out;
        out << q_.R(), t_, (T)0., (T)0., (T)0., (T)1.;
        return out;
    }

    /**
     * @brief Obtain the inverse transform \f$\boldsymbol{T}_A^B\rightarrow \boldsymbol{T}_B^A\f$.
     */
    SE3 inverse() const
    {
        SE3    x;
        SO3<T> q_inv = q_.inverse();
        return SE3::fromVecAndQuat(-(q_inv * t_), q_inv);
    }

    /**
     * @brief Invert the current transform \f$\boldsymbol{T}_A^B\rightarrow \boldsymbol{T}_B^A\f$.
     */
    SE3& invert()
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
    SE3<Tout> otimes(const SE3<T2>& x) const
    {
        SE3<Tout> xout;
        xout.t() = t_ + q_ * x.t_;
        xout.q() = q_ * x.q_;
        return xout;
    }

    /**
     * @brief Implementation of tangent space group perturbations: \f$\boldsymbol{T}_A^B\oplus \boldsymbol{t}_B^{B'}
     * \rightarrow \boldsymbol{T}_A^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    SE3<Tout> oplus(const Matrix<T2, 6, 1>& delta) const
    {
        return otimes<Tout, T2>(SE3<T2>::Exp(delta));
    }

    /**
     * @brief Implementation of group subtraction: \f$\boldsymbol{T}_A^B\ominus \boldsymbol{T}_A^{B'} \rightarrow
     * \boldsymbol{t}_B^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 6, 1> ominus(const SE3<T2>& x) const
    {
        return SE3<Tout>::Log(x.inverse().template otimes<Tout>(*this));
    }

    /**
     * @brief Copy constructor.
     */
    SE3& operator=(const SE3& x)
    {
        arr_ = x.elements();
        return *this;
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SE3 operator*(const SE3<T2>& x) const
    {
        return otimes(x);
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SE3& operator*=(const SE3<T2>& x)
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
    SE3& operator*=(const double& s)
    {
        arr_ = SE3::Exp(s * SE3::Log(*this)).elements();
        return *this;
    }

    /**
     * @brief Scale a transform by a scalar.
     *
     * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a transform.
     */
    SE3& operator/=(const double& s)
    {
        arr_ = SE3::Exp(SE3::Log(*this) / s).elements();
        return *this;
    }

    /**
     * @brief Scale a transform by a scalar.
     *
     * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a transform.
     */
    SE3 operator/(const double& s) const
    {
        SE3 qs;
        qs.arr_ = SE3::Exp(SE3::Log(*this) / s).elements();
        return qs;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{T}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 3, 1> operator*(const Matrix<T2, 3, 1>& v) const
    {
        return q_ * v + t_;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{T}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    Vec3T operator*(const Vec3T& v) const
    {
        return q_ * v + t_;
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SE3 operator+(const Vec6T& v) const
    {
        return oplus(v);
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SE3& operator+=(const Vec6T& v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }

    /**
     * @brief Invocation of ominus via subtraction.
     */
    template<typename T2>
    Vec6T operator-(const SE3<T2>& x) const
    {
        return ominus(x);
    }

    /**
     * @brief Hat operator implementation, which coverts the tangent-space vector representation to the corresponding
     * Lie algebra: \f$\mathbb{R}^6\rightarrow \mathfrak{se}(3)\f$.
     */
    static Mat4T hat(const Vec6T& omega)
    {
        Mat4T Omega;
        Omega << SO3<T>::hat(omega.template block<3, 1>(3, 0)), omega.template block<3, 1>(0, 0), (T)0., (T)0., (T)0.,
            (T)0.;
        return Omega;
    }

    /**
     * @brief Vee operator implementation, which coverts the Lie algebra representation to a tangent-space vector
     * representation: \f$\mathfrak{se}(3) \rightarrow \mathbb{R}^6\f$.
     */
    static Vec6T vee(const Mat4T& Omega)
    {
        Vec6T omega;
        omega << Omega.template block<3, 1>(0, 3), SO3<T>::vee(Omega.template block<3, 3>(0, 0));
        return omega;
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SE(3) \rightarrow \mathfrak{se}(3)\f$.
     */
    static Mat4T log(const SE3& x)
    {
        return SE3::hat(SE3::Log(x));
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SE(3) \rightarrow \mathbb{R}^6\f$.
     */
    static Vec6T Log(const SE3& x)
    {
        Vec3T w  = SO3<T>::Log(x.q_);
        T     th = w.norm();
        Mat3T W  = SO3<T>::hat(w);

        Mat3T leftJacobianInverse;
        if (th > (T)1e-4)
        {
            T a                 = sin(th) / th;
            T b                 = ((T)1. - cos(th)) / (th * th);
            T c                 = ((T)1. - a) / (th * th);
            T e                 = (b - 2. * c) / (2. * a);
            leftJacobianInverse = Mat3T::Identity() - 0.5 * W + e * (W * W);
        }
        else
        {
            leftJacobianInverse = Mat3T::Identity();
        }

        Vec6T out;
        out << leftJacobianInverse * x.t_, w;

        return out;
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathfrak{se}(3) \rightarrow SE(3)\f$.
     */
    static SE3 exp(const Mat4T& Omega)
    {
        return SE3::Exp(SE3::vee(Omega));
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathbb{R}^6 \rightarrow SE(3)\f$.
     */
    static SE3 Exp(const Vec6T& omega)
    {
        Vec3T  rho = omega.template block<3, 1>(0, 0);
        Vec3T  w   = omega.template block<3, 1>(3, 0);
        SO3<T> q   = SO3<T>::Exp(w);
        Mat3T  W   = SO3<T>::hat(w);
        T      th  = w.norm();

        Mat3T leftJacobian;
        if (th > (T)1e-4)
        {
            T a          = sin(th) / th;
            T b          = ((T)1. - cos(th)) / (th * th);
            T c          = ((T)1. - a) / (th * th);
            leftJacobian = a * Mat3T::Identity() + b * W + c * (w * w.transpose());
        }
        else
        {
            leftJacobian = Mat3T::Identity();
        }

        return SE3::fromVecAndQuat(leftJacobian * rho, q);
    }

    /**
     * @brief Cast the underlying numeric type.
     */
    template<typename T2>
    SE3<T2> cast() const
    {
        SE3<T2> x;
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
SE3<T> operator*(const double& l, const SE3<T>& r)
{
    SE3<T> lr;
    lr.arr_ = SE3<T>::Exp(l * SE3<T>::Log(r)).elements();
    return lr;
}

/**
 * @brief Scale a transform by a scalar.
 *
 * Under the hood, this converts the transform into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a transform.
 */
template<typename T>
SE3<T> operator*(const SE3<T>& l, const double& r)
{
    SE3<T> lr;
    lr.arr_ = SE3<T>::Exp(r * SE3<T>::Log(l)).elements();
    return lr;
}

/**
 * @brief Render the transform in a stream.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SE3<T>& x)
{
    os << "SE(3): [ " << x.t_.x() << "i, " << x.t_.y() << "j, " << x.t_.z() << "k ] [ " << x.q_.w() << ", " << x.q_.x()
       << "i, " << x.q_.y() << "j, " << x.q_.z() << "k ]";
    return os;
}

typedef SE3<double> SE3d;
