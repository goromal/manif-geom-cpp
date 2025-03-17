#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <math.h>

using namespace Eigen;

/**
 * @brief Class representing a member of the \f$SO(3)\f$ manifold, or a 3D rotation.
 */
template<typename T>
class SO3
{
private:
    typedef Matrix<T, 3, 1> Vec3T;
    typedef Matrix<T, 4, 1> Vec4T;
    typedef Matrix<T, 3, 3> Mat3T;
    typedef Matrix<T, 4, 4> Mat4T;

    T buf_[4];

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /**
     * @brief Memory-mapped array representing all rotation fields in \f$\boldsymbol{q}\f$.
     */
    Map<Vec4T> arr_;

    /**
     * @brief Obtain a random rotation.
     * @return A random rotation \f$\boldsymbol{q}_B^W\in SO(3)\f$.
     *
     *  The rotation \f$\mathbf{q}_B^W\f$ will be normalized.
     */
    static SO3 random()
    {
        SO3 q;
        q.arr_.setRandom();
        q.normalize();
        return q;
    }

    /**
     * @brief Obtain an identity \f$SO(3)\f$ rotation.
     */
    static SO3 identity()
    {
        SO3 q;
        q.arr_ << 1., 0., 0., 0.;
        return q;
    }

    /**
     * @brief Obtain a rotation full of NaNs.
     */
    static SO3 nans()
    {
        SO3 x;
        x.arr_.setConstant(std::numeric_limits<T>::quiet_NaN());
        return x;
    }

    /**
     * @brief Convert an angle (in radians) with accompanying axis (defined by a vector) into a rotation.
     */
    static SO3 fromAxisAngle(const Vec3T& axis, const T angle)
    {
        T     th2 = angle / 2.0;
        Vec3T xyz = sin(th2) * axis / axis.norm();
        SO3   q;
        q.arr_ << cos(th2), xyz(0), xyz(1), xyz(2);
        q.normalize();
        return q;
    }

    /**
     * @brief Construct a rotation from yaw-pitch-roll successive-axes Euler angles.
     */
    static SO3 fromEuler(const T roll, const T pitch, const T yaw)
    {
        SO3 q_roll  = SO3::fromAxisAngle((Vec3T() << (T)1., (T)0., (T)0.).finished(), roll);
        SO3 q_pitch = SO3::fromAxisAngle((Vec3T() << (T)0., (T)1., (T)0.).finished(), pitch);
        SO3 q_yaw   = SO3::fromAxisAngle((Vec3T() << (T)0., (T)0., (T)1.).finished(), yaw);
        SO3 q_euler = q_yaw * q_pitch * q_roll;
        q_euler.normalize();
        return q_euler;
    }

    /**
     * @brief Convert a rotation matrix \f$\boldsymbol{R}\in\mathbb{R}^{3\times 3}\f$ to a \f$SO(3)\f$ rotation.
     */
    static SO3 fromR(const Mat3T& m)
    {
        T s, qw, qx, qy, qz;

        T R11 = m(0, 0);
        T R12 = m(0, 1);
        T R13 = m(0, 2);
        T R21 = m(1, 0);
        T R22 = m(1, 1);
        T R23 = m(1, 2);
        T R31 = m(2, 0);
        T R32 = m(2, 1);
        T R33 = m(2, 2);

        if (R11 + R22 + R33 > (T)0.0)
        {
            s  = 2. * sqrt(1. + R11 + R22 + R33);
            qw = 0.25 * s;
            qx = (R32 - R23) / s;
            qy = (R13 - R31) / s;
            qz = (R21 - R12) / s;
        }
        else if (R11 > R22 && R11 > R33)
        {
            s  = 2. * sqrt(1. + R11 - R22 - R33);
            qw = (R32 - R23) / s;
            qx = 0.25 * s;
            qy = (R21 + R12) / s;
            qz = (R31 + R13) / s;
        }
        else if (R22 > R33)
        {
            s  = 2. * sqrt(1. + R22 - R11 - R33);
            qw = (R13 - R31) / s;
            qx = (R21 + R12) / s;
            qy = 0.25 * s;
            qz = (R32 + R23) / s;
        }
        else
        {
            s  = 2. * sqrt(1. + R33 - R11 - R22);
            qw = (R21 - R12) / s;
            qx = (R31 + R13) / s;
            qy = (R32 + R23) / s;
            qz = 0.25 * s;
        }

        SO3 q = SO3::fromQuat(qw, qx, qy, qz);
        q.normalize();

        return q;
    }

    /**
     * @brief Given two unit vectors \f$\boldsymbol{u},\boldsymbol{v}\in \mathbb{R}^3\f$, returns the rotation that
     * rotates \f$\boldsymbol{u}\rightarrow\boldsymbol{v}\f$.
     *
     * Naturally, there is roll ambiguity in such a transform, so watch out.
     */
    static SO3 fromTwoUnitVectors(const Vec3T& u, const Vec3T& v)
    {
        SO3 q;
        T   d = u.dot(v);

        if (d < (T)0.99999999 && d > (T)-0.99999999)
        {
            T     invs = 1. / sqrt((2. * (1. + d)));
            Vec3T xyz  = u.cross(v * invs);
            q          = SO3::fromQuat(0.5 / invs, xyz(0), xyz(1), xyz(2));
            q.normalize();
        }
        else if (d < (T)-0.99999999)
        {
            // There are an infinite number of solutions here, choose one.
            // This choice works better for vector comparisons with only
            // nonzero x components.
            q = SO3::fromQuat((T)0., (T)0., (T)1., (T)0.);
        }
        else
        {
            q = SO3::identity();
        }

        return q;
    }

    /**
     * @brief Construct a rotation from the individual fields.
     * @param qw Real component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qx First imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qy Second imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     * @param qz Third imaginary component of the rotation \f$\boldsymbol{q}\in SO(3)\f$.
     */
    static SO3 fromQuat(const T qw, const T qx, const T qy, const T qz)
    {
        SO3 q;
        q.arr_ << qw, qx, qy, qz;
        return q;
    }

    /**
     * @brief Construct a rotation from a rotation fields vector.
     * @param qvec The rotation represented as an array \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    static SO3 fromQuat(const Vec4T& qvec)
    {
        SO3 q;
        q.arr_ << qvec(0), qvec(1), qvec(2), qvec(3);
        return q;
    }

    /**
     * @brief Construct a rotation from an Eigen quaternion.
     * @param quat The Eigen quaternion.
     */
    static SO3 fromQuat(const Quaternion<T> quat)
    {
        SO3 q;
        q.arr_ << quat.w(), quat.x(), quat.y(), quat.z();
        return q;
    }

    /**
     * @brief Create a rotation (with garbage data).
     */
    SO3() : arr_(buf_) {}

    /**
     * @brief Create a transform from an array representing all rotation fields in \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    SO3(const Ref<const Vec4T>& arr) : arr_(buf_)
    {
        arr_ = arr;
    }

    /**
     * @brief Copy constructor from another rotation.
     */
    SO3(const SO3& q) : arr_(buf_)
    {
        arr_ = q.arr_;
    }

    /**
     * @brief Create a rotation from a pointer array representing all rotation fields in
     * \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    SO3(const T* data) : arr_(const_cast<T*>(data)) {}

    /**
     * @brief Access a field from \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
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
     * @brief Access the first imaginary element of the rotation.
     */
    inline const T& x() const
    {
        return arr_(1);
    }

    /**
     * @brief Access the second imaginary element of the rotation.
     */
    inline const T& y() const
    {
        return arr_(2);
    }

    /**
     * @brief Access the third imaginary element of the rotation.
     */
    inline const T& z() const
    {
        return arr_(3);
    }

    /**
     * @brief Access the real element of the rotation.
     */
    inline T& w()
    {
        return arr_(0);
    }

    /**
     * @brief Access the first imaginary element of the rotation.
     */
    inline T& x()
    {
        return arr_(1);
    }

    /**
     * @brief Access the second imaginary element of the rotation.
     */
    inline T& y()
    {
        return arr_(2);
    }

    /**
     * @brief Access the third imaginary element of the rotation.
     */
    inline T& z()
    {
        return arr_(3);
    }

    /**
     * @brief Access all elements of \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    inline const Vec4T elements() const
    {
        return arr_;
    }

    /**
     * @brief Access all elements of \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    inline Vec4T array() const
    {
        return arr_;
    }

    /**
     * @brief Access pointer to all elements of \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    inline T* data()
    {
        return arr_.data();
    }

    /**
     * @brief Access pointer to all elements of \f$\boldsymbol{q}\in\mathbb{R}^4\f$.
     */
    inline const T* data() const
    {
        return arr_.data();
    }

    /**
     * @brief Get a deep copy of the current rotation.
     */
    SO3 copy() const
    {
        SO3 tmp;
        tmp.arr_ = arr_;
        return tmp;
    }

    /**
     * @brief Normalize the elements of \f$\boldsymbol{q}\f$ to make it a valid rotation.
     */
    void normalize()
    {
        arr_ /= arr_.norm();
        if (arr_(0) < (T)0.0)
            arr_ *= (T)-1.0;
    }

    /**
     * @brief Obtain a normalized copy of the current rotation.
     */
    SO3 normalized()
    {
        SO3 tmp = copy();
        tmp.normalize();
        return tmp;
    }

    /**
     * @brief Convert the rotation to matrix representation \f$\boldsymbol{R}\in \mathbb{R}^{3\times 3}\f$.
     */
    Mat3T R() const
    {
        T     wx = w() * x();
        T     wy = w() * y();
        T     wz = w() * z();
        T     xx = x() * x();
        T     xy = x() * y();
        T     xz = x() * z();
        T     yy = y() * y();
        T     yz = y() * z();
        T     zz = z() * z();
        Mat3T out;
        out << 1. - 2. * yy - 2. * zz, 2. * xy - 2. * wz, 2. * xz + 2. * wy, 2. * xy + 2. * wz, 1. - 2. * xx - 2. * zz,
            2. * yz - 2. * wx, 2. * xz - 2. * wy, 2. * yz + 2. * wx, 1. - 2. * xx - 2. * yy;
        return out;
    }

    /**
     * @brief Obtain the inverse rotation \f$\boldsymbol{R}_A^B\rightarrow \boldsymbol{R}_B^A\f$.
     */
    SO3 inverse() const
    {
        SO3 q;
        q.arr_(0) = arr_(0);
        q.arr_(1) = -arr_(1);
        q.arr_(2) = -arr_(2);
        q.arr_(3) = -arr_(3);
        return q;
    }

    /**
     * @brief Invert the current rotation \f$\boldsymbol{R}_A^B\rightarrow \boldsymbol{R}_B^A\f$.
     */
    SO3& invert()
    {
        arr_.template block<3, 1>(1, 0) *= (T)-1.0;
        return *this;
    }

    /**
     * @brief Obtain the roll component of the yaw-pitch-roll successive-axes Euler angle representation.
     */
    T roll() const
    {
        return atan2(T(2.0) * (w() * x() + y() * z()), T(1.0) - T(2.0) * (x() * x() + y() * y()));
    }

    /**
     * @brief Obtain the pitch component of the yaw-pitch-roll successive-axes Euler angle representation.
     */
    T pitch() const
    {
        const T val = T(2.0) * (w() * y() - x() * z());

        // hold at 90 degrees if invalid
        if (fabs(val) > T(1.0))
            return copysign(T(1.0), val) * T(M_PI) / T(2.0);
        else
            return asin(val);
    }

    /**
     * @brief Obtain the yaw component of the yaw-pitch-roll successive-axes Euler angle representation.
     */
    T yaw() const
    {
        return atan2(T(2.0) * (w() * z() + x() * y()), T(1.0) - T(2.0) * (y() * y() + z() * z()));
    }

    /**
     * @brief Obtain the Euler angle representation as a 3-vector.
     */
    Vec3T toEuler() const
    {
        Vec3T out;
        out << roll(), pitch(), yaw();
        return out;
    }

    /**
     * @brief Obtain an equivalent 4-by-4 matrix representation of \f$\boldsymbol{q}\rightarrow [\boldsymbol{q}]_L\f$
     * such that \f$\boldsymbol{q}\otimes \boldsymbol{q}_{\text{other}}=[\boldsymbol{q}]_L
     * \boldsymbol{q}_{\text{other}}\f$.
     */
    Mat4T qMatLeft() const
    {
        Mat4T qL;
        qL << w(), -x(), -y(), -z(), x(), w(), -z(), y(), y(), z(), w(), -x(), z(), -y(), x(), w();
        return qL;
    }

    /**
     * @brief Implementation of group composition: \f$\boldsymbol{q}_B^C \otimes \boldsymbol{q}_A^B\rightarrow
     * \boldsymbol{q}_A^C\f$.
     */
    template<typename Tout = T, typename T2>
    SO3<Tout> otimes(const SO3<T2>& q) const
    {
        SO3<Tout> qout;
        qout.arr_ << w() * q.w() - x() * q.x() - y() * q.y() - z() * q.z(),
            w() * q.x() + x() * q.w() + y() * q.z() - z() * q.y(),
            w() * q.y() - x() * q.z() + y() * q.w() + z() * q.x(),
            w() * q.z() + x() * q.y() - y() * q.x() + z() * q.w();
        return qout;
    }

    /**
     * @brief Implementation of tangent space group perturbations: \f$\boldsymbol{q}_A^B\oplus
     * \boldsymbol{\theta}_B^{B'} \rightarrow \boldsymbol{q}_A^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    SO3<Tout> oplus(const Matrix<T2, 3, 1>& delta) const
    {
        return otimes<Tout, T2>(SO3<T2>::Exp(delta));
    }

    /**
     * @brief Implementation of group subtraction: \f$\boldsymbol{q}_A^B\ominus \boldsymbol{q}_A^{B'} \rightarrow
     * \boldsymbol{\theta}_B^{B'}\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 3, 1> ominus(const SO3<T2>& q) const
    {
        SO3<Tout> dq = q.inverse().template otimes<Tout>(*this);
        if (dq.w() < 0.0)
        {
            dq.arr_ *= (Tout)-1.0;
        }
        return SO3<Tout>::Log(dq);
    }

    /**
     * @brief Copy constructor.
     */
    SO3& operator=(const SO3& q)
    {
        arr_ = q.elements();
        return *this;
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SO3 operator*(const SO3<T2>& q) const
    {
        return otimes(q);
    }

    /**
     * @brief Invocation of otimes via multiplication.
     */
    template<typename T2>
    SO3& operator*=(const SO3<T2>& q)
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
    SO3& operator*=(const double& s)
    {
        arr_ = SO3::Exp(s * SO3::Log(*this)).elements();
        return *this;
    }

    /**
     * @brief Scale a rotation by a scalar.
     *
     * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a rotation.
     */
    SO3& operator/=(const double& s)
    {
        arr_ = SO3::Exp(SO3::Log(*this) / s).elements();
        return *this;
    }

    /**
     * @brief Scale a rotation by a scalar.
     *
     * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
     * scaled vector back to a rotation.
     */
    SO3 operator/(const double& s) const
    {
        SO3 qs;
        qs.arr_ = SO3::Exp(SO3::Log(*this) / s).elements();
        return qs;
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{q}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    template<typename Tout = T, typename T2>
    Matrix<Tout, 3, 1> operator*(const Matrix<T2, 3, 1>& v) const
    {
        Vec3T              qv = arr_.template block<3, 1>(1, 0);
        Matrix<Tout, 3, 1> t  = (Tout)2.0 * v.cross(qv);
        return v - w() * t + t.cross(qv);
    }

    /**
     * @brief Transform a vector via multiplication: \f$\boldsymbol{q}_A^B\boldsymbol{t}^A \rightarrow
     * \boldsymbol{t}^B\f$.
     */
    Vec3T operator*(const Vec3T& v) const
    {
        Vec3T qv = arr_.template block<3, 1>(1, 0);
        Vec3T t  = (T)2.0 * v.cross(qv);
        return v - w() * t + t.cross(qv);
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SO3 operator+(const Vec3T& v) const
    {
        return oplus(v);
    }

    /**
     * @brief Invocation of oplus via addition.
     */
    SO3& operator+=(const Vec3T& v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }

    /**
     * @brief Invocation of ominus via subtraction.
     */
    template<typename T2>
    Vec3T operator-(const SO3<T2>& q) const
    {
        return ominus(q);
    }

    /**
     * @brief Hat operator implementation, which coverts the tangent-space vector representation to the corresponding
     * Lie algebra: \f$\mathbb{R}^3\rightarrow \mathfrak{so}(3)\f$.
     */
    static Mat3T hat(const Vec3T& omega)
    {
        Mat3T Omega;
        Omega << (T)0., -omega.z(), omega.y(), omega.z(), (T)0., -omega.x(), -omega.y(), omega.x(), (T)0.;
        return Omega;
    }

    /**
     * @brief Vee operator implementation, which coverts the Lie algebra representation to a tangent-space vector
     * representation: \f$\mathfrak{so}(3) \rightarrow \mathbb{R}^3\f$.
     */
    static Vec3T vee(const Mat3T& Omega)
    {
        Vec3T omega;
        omega << Omega(2, 1), Omega(0, 2), Omega(1, 0);
        return omega;
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SO(3) \rightarrow \mathfrak{so}(3)\f$.
     */
    static Mat3T log(const SO3& q)
    {
        return hat(SO3::Log(q));
    }

    /**
     * @brief Logarithmic chart map implementation: \f$SO(3) \rightarrow \mathbb{R}^3\f$.
     */
    static Vec3T Log(const SO3& q)
    {
        Vec3T qv = q.elements().template block<3, 1>(1, 0);
        T     qw = q.w();

        T n = qv.norm();
        if (n > (T)1e-4)
            return 2.0 * qv * atan2(n, qw) / n;
        else
            return qv;
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathfrak{so}(3) \rightarrow SO(3)\f$.
     */
    static SO3 exp(const Mat3T& Omega)
    {
        return SO3::Exp(vee(Omega));
    }

    /**
     * @brief Exponential chart map implementation: \f$\mathbb{R}^3 \rightarrow SO(3)\f$.
     */
    static SO3 Exp(const Vec3T& omega)
    {
        T th = omega.norm();

        SO3 q;
        if (th > (T)1e-4)
        {
            Vec3T u                           = omega / th;
            q.arr_(0)                         = cos(th / 2.0);
            q.arr_.template block<3, 1>(1, 0) = sin(th / 2.0) * u;
        }
        else
        {
            q.arr_(0)                         = (T)1.0;
            q.arr_.template block<3, 1>(1, 0) = omega / 2.0;
            q.normalize();
        }

        return q;
    }

    /**
     * @brief Cast the underlying numeric type.
     */
    template<typename T2>
    SO3<T2> cast() const
    {
        SO3<T2> q;
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
SO3<T> operator*(const double& l, const SO3<T>& r)
{
    SO3<T> lr;
    lr.arr_ = SO3<T>::Exp(l * SO3<T>::Log(r)).elements();
    return lr;
}

/**
 * @brief Scale a rotation by a scalar.
 *
 * Under the hood, this converts the rotation into a tangent-space vector, scales the vector, then converts the
 * scaled vector back to a rotation.
 */
template<typename T>
SO3<T> operator*(const SO3<T>& l, const double& r)
{
    SO3<T> lr;
    lr.arr_ = SO3<T>::Exp(r * SO3<T>::Log(l)).elements();
    return lr;
}

/**
 * @brief Render the rotation in a stream.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const SO3<T>& q)
{
    os << "SO(3): [ " << q.w() << ", " << q.x() << "i, " << q.y() << "j, " << q.z() << "k ]";
    return os;
}

typedef SO3<double> SO3d;
