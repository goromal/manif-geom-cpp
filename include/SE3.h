#pragma once

#include "SO3.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <math.h>

using namespace Eigen;

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
    Map<Vec7T> arr_;
    Map<Vec3T> t_;
    SO3<T>     q_;

    static SE3 random()
    {
        SE3 x;
        x.arr_.setRandom();
        x.q_.normalize();
        return x;
    }

    static SE3 identity()
    {
        SE3 x;
        x.arr_.setZero();
        x.arr_(3) = (T)1.0;
        return x;
    }

    static SE3 nans()
    {
        SE3 x;
        x.arr_.setConstant(std::numeric_limits<T>::quiet_NaN());
        return x;
    }

    static SE3 fromH(const Mat4T& m)
    {
        return SE3::fromVecAndQuat(m.template block<3, 1>(0, 3), SO3<T>::fromR(m.template block<3, 3>(0, 0)));
    }

    static SE3 fromVecAndQuat(const T tx, const T ty, const T tz, const T qw, const T qx, const T qy, const T qz)
    {
        SE3 x;
        x.arr_ << tx, ty, tz, qw, qx, qy, qz;
        return x;
    }

    static SE3 fromVecAndQuat(const Vec3T& tvec, const Vec4T& qvec)
    {
        SE3 x;
        x.arr_ << tvec(0), tvec(1), tvec(2), qvec(0), qvec(1), qvec(2), qvec(3);
        return x;
    }

    static SE3 fromVecAndQuat(const Vec3T& t, const SO3<T>& q)
    {
        SE3 x;
        x.arr_ << t(0), t(1), t(2), q.w(), q.x(), q.y(), q.z();
        return x;
    }

    static SE3 fromVecAndQuat(const Vec3T& tvec, const Quaternion<T> quat)
    {
        SE3 x;
        x.arr_ << tvec(0), tvec(1), tvec(2), quat.w(), quat.x(), quat.y(), quat.z();
        return x;
    }

    SE3() : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3) {}

    SE3(const Ref<const Vec7T>& arr) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3)
    {
        arr_ = arr;
    }

    SE3(const SE3& x) : arr_(buf_), t_(arr_.data()), q_(arr_.data() + 3)
    {
        arr_ = x.arr_;
    }

    SE3(const T* data) : arr_(const_cast<T*>(data)), t_(arr_.data()), q_(arr_.data() + 3) {}

    inline T& operator[](int i)
    {
        return arr_[i];
    }
    inline const Map<Vec3T>& t() const
    {
        return t_;
    }
    inline const SO3<T>& q() const
    {
        return q_;
    }
    inline Map<Vec3T>& t()
    {
        return t_;
    }
    inline SO3<T>& q()
    {
        return q_;
    }
    inline const Vec7T elements() const
    {
        return arr_;
    }
    inline Vec7T array() const
    {
        return arr_;
    }
    inline T* data()
    {
        return arr_.data();
    }
    inline const T* data() const
    {
        return arr_.data();
    }

    SE3 copy() const
    {
        SE3 tmp;
        tmp.arr_ = arr_;
        return tmp;
    }

    Mat4T H() const
    {
        Mat4T out;
        out << q_.R(), t_, (T)0., (T)0., (T)0., (T)1.;
        return out;
    }

    SE3 inverse() const
    {
        SE3    x;
        SO3<T> q_inv = q_.inverse();
        return SE3::fromVecAndQuat(-(q_inv * t_), q_inv);
    }

    SE3& invert()
    {
        q().invert();
        t() = -(q() * t());
        return *this;
    }

    template<typename Tout = T, typename T2>
    SE3<Tout> otimes(const SE3<T2>& x) const
    {
        SE3<Tout> xout;
        xout.t() = t_ + q_ * x.t_;
        xout.q() = q_ * x.q_;
        return xout;
    }

    template<typename Tout = T, typename T2>
    SE3<Tout> oplus(const Matrix<T2, 6, 1>& delta) const
    {
        return otimes<Tout, T2>(SE3<T2>::Exp(delta));
    }

    template<typename Tout = T, typename T2>
    Matrix<Tout, 6, 1> ominus(const SE3<T2>& x) const
    {
        return SE3<Tout>::Log(x.inverse().template otimes<Tout>(*this));
    }

    SE3& operator=(const SE3& x)
    {
        arr_ = x.elements();
        return *this;
    }

    template<typename T2>
    SE3 operator*(const SE3<T2>& x) const
    {
        return otimes(x);
    }

    template<typename T2>
    SE3& operator*=(const SE3<T2>& x)
    {
        arr_ = otimes(x).elements();
        return *this;
    }

    SE3& operator*=(const double& s)
    {
        arr_ = SE3::Exp(s * SE3::Log(*this)).elements();
        return *this;
    }

    SE3& operator/=(const double& s)
    {
        arr_ = SE3::Exp(SE3::Log(*this) / s).elements();
        return *this;
    }

    SE3 operator/(const double& s) const
    {
        SE3 qs;
        qs.arr_ = SE3::Exp(SE3::Log(*this) / s).elements();
        return qs;
    }

    template<typename Tout = T, typename T2>
    Matrix<Tout, 3, 1> operator*(const Matrix<T2, 3, 1>& v) const
    {
        return q_ * v + t_;
    }

    Vec3T operator*(const Vec3T& v) const
    {
        return q_ * v + t_;
    }

    SE3 operator+(const Vec6T& v) const
    {
        return oplus(v);
    }

    SE3& operator+=(const Vec6T& v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }

    template<typename T2>
    Vec6T operator-(const SE3<T2>& x) const
    {
        return ominus(x);
    }

    static Mat4T hat(const Vec6T& omega)
    {
        Mat4T Omega;
        Omega << SO3<T>::hat(omega.template block<3, 1>(3, 0)), omega.template block<3, 1>(0, 0), (T)0., (T)0., (T)0.,
            (T)0.;
        return Omega;
    }

    static Vec6T vee(const Mat4T& Omega)
    {
        Vec6T omega;
        omega << Omega.template block<3, 1>(0, 3), SO3<T>::vee(Omega.template block<3, 3>(0, 0));
        return omega;
    }

    static Mat4T log(const SE3& x)
    {
        return SE3::hat(SE3::Log(x));
    }

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

    static SE3 exp(const Mat4T& Omega)
    {
        return SE3::Exp(SE3::vee(Omega));
    }

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

    template<typename T2>
    SE3<T2> cast() const
    {
        SE3<T2> x;
        x.arr_ = arr_.template cast<T2>();
        return x;
    }
};

template<typename T>
SE3<T> operator*(const double& l, const SE3<T>& r)
{
    SE3<T> lr;
    lr.arr_ = SE3<T>::Exp(l * SE3<T>::Log(r)).elements();
    return lr;
}

template<typename T>
SE3<T> operator*(const SE3<T>& l, const double& r)
{
    SE3<T> lr;
    lr.arr_ = SE3<T>::Exp(r * SE3<T>::Log(l)).elements();
    return lr;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SE3<T>& x)
{
    os << "SE(3): [ " << x.t_.x() << "i, " << x.t_.y() << "j, " << x.t_.z() << "k ] [ " << x.q_.w() << ", " << x.q_.x()
       << "i, " << x.q_.y() << "j, " << x.q_.z() << "k ]";
    return os;
}

typedef SE3<double> SE3d;
