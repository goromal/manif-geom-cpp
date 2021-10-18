#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>
#include <iostream>

using namespace Eigen;

template<typename T>
class SO3
{

private:
  typedef Matrix<T,3,1> Vec3T;
  typedef Matrix<T,4,1> Vec4T;
  typedef Matrix<T,3,3> Mat3T;
  typedef Matrix<T,4,4> Mat4T;
  
  static const Vec3T e1 = [] {
    Vec3T v;
    v << (T)1, (T)0, (T)0;
    return v;
  }();
  
  static const Vec3T e2 = [] {
    Vec3T v;
    v << (T)0, (T)1, (T)0;
    return v;
  }();
  
  static const Vec3T e3 = [] {
    Vec3T v;
    v << (T)0, (T)0, (T)1;
    return v;
  }();
  
  T buf_[4];

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Map<Vec4T> arr_;
  
  static SO3 random()
  {
    SO3 q;
    q.arr_.setRandom();
    q.normalize();
    return q;
  }
  
  static SO3 identity()
  {
    SO3 q;
    q.arr_ << 1.0, 0, 0, 0;
    return q;
  }
  
  static SO3 fromAxisAngle(const Vec3T &axis, const T angle)
  {
    T th2 = angle / 2.0;
    SO3 q;
    q.arr_ << cos(th2), sin(th2)*axis/axis.norm();
    q.normalize();
    return q;
  }
  
  static SO3 fromEuler(const T roll, const T pitch, const T yaw)
  {
    SO3 q_roll  = SO3.fromAxisAngle(e1, roll);
    SO3 q_pitch = SO3.fromAxisAngle(e2, pitch);
    SO3 q_yaw   = SO3.fromAxisAngle(e3, yaw);
    return q_yaw * q_pitch * q_roll;
  }
  
  static SO3 fromR(const Mat3T &m)
  {
    T s, qw, qx, qy, qz;
  
    T R11 = m(0,0);
    T R12 = m(0,1);
    T R13 = m(0,2);
    T R21 = m(1,0);
    T R22 = m(1,1);
    T R23 = m(1,2);
    T R31 = m(2,0);
    T R32 = m(2,1);
    T R33 = m(2,2);
    
    if (R11 + R22 + R33 > 0)
    {
      s  = 2. * sqrt(1. + R11 + R22 + R33);
      qw = 0.25 * s;
      qx = (R32-R23) / s;
      qy = (R13-R31) / s;
      qz = (R21-R12) / s;
    }
    else if (R11 > R22 && R11 > R33)
    {
      s  = 2. * sqrt(1. + R11 - R22 - R33);
      qw = (R32-R23) / s;
      qx = 0.25 * s;
      qy = (R21+R12) / s;
      qz = (R31+R13) / s;
    }
    else if (R22 > R33)
    {
      s  = 2. * sqrt(1. + R22 - R11 - R33);
      qw = (R13-R31) / s;
      qx = (R21+R12) / s;
      qy = 0.25 * s;
      qz = (R32+R23) / s;
    }
    else
    {
      s  = 2. * sqrt(1. + R33 - R11 - R22);
      qw = (R21-R12) / s;
      qx = (R31+R13) / s;
      qy = (R32+R23) / s;
      qz = 0.25 * s;
    }
        
    SO3 q = SO3.fromQuat(qw, qx, qy, qz);
    q.normalize();

    return q;
  }
  
  static SO3 = fromTwoUnitVectors(const Vec3T &u, const Vec3T &v)
  {
    SO3 q;
    T d = u.dot(v);
    
    if (d < 0.99999999 && d > -0.99999999)
    {
      T invs = 1. / sqrt((2.*(1.+d)));
      Vec3 xyz = u.cross(v*invs);
      q = SO3.fromQuat(0.5/invs, xyz(0), xyz(1), xyz(2));
      q.normalize();
    }
    else if (d < -0.99999999)
    {
      // There are an infinite number of solutions here, choose one.
      // This choice works better for vector comparisons with only 
      // nonzero x components.
      q = SO3.fromQuat((T)0, (T)0, (T)1, (T)0);
    }
    else
    {
      q = SO3.identity();
    }
    
    return q;
  }
  
  static SO3 fromQuat(const T qw, const T qx, const T qy, const T qz)
  {
    SO3 q;
    q.arr_ << qw, qx, qy, qz;
    return q;
  }
  
  static SO3 fromQuat(const Vec4T &qvec)
  {
    SO3 q;
    q.arr_ << qvec(0), qvec(1), qvec(2), qvec(3);
    return q;
  }
  
  static SO3 fromQuat(const Quaternion<T> quat)
  {
    SO3 q;
    q.arr_ << quat.w(), quat.x(), quat.y(), quat.z();
    return q;
  }
  
  SO3() :
    arr_(buf_)
  {}

  SO3(const Ref<const Vec4T> arr) :
    arr_(const_cast<T*>(arr.data()))
  {}

  SO3(const SO3& q)
    : arr_(buf_)
  {
    arr_ = q.arr_;
  }

  SO3(const T* data) :
    arr_(const_cast<T*>(data))
  {}

  inline T* data() { return arr_.data(); 
  inline const T* data() const { return arr_.data();}
  inline T& operator[] (int i) {return arr_[i];}
  inline T w() const { return arr_(0); }
  inline T x() const { return arr_(1); }
  inline T y() const { return arr_(2); }
  inline T z() const { return arr_(3); }
  inline void setW(T w) { arr_(0) = w; }
  inline void setX(T x) { arr_(1) = x; }
  inline void setY(T y) { arr_(2) = y; }
  inline void setZ(T z) { arr_(3) = z; }
  inline const Vec4T array() const { return arr_;}
  
  void normalize()
  {
    arr_ /= arr_.norm();
    if (arr_(0) < 0) arr_ *= (T)-1.0;
  }
  
  Mat3T R() const
  {
    T wx = w()*x();
    T wy = w()*y();
    T wz = w()*z();
    T xx = x()*x();
    T xy = x()*y();
    T xz = x()*z();
    T yy = y()*y();
    T yz = y()*z();
    T zz = z()*z();
    Matrix3T out;
    out << 1. - 2.*yy - 2.*zz, 2.*xy - 2.*wz,      2.*xz + 2.*wy,
           2.*xy + 2.*wz,      1. - 2.*xx - 2.*zz, 2.*yz - 2.*wx,
           2.*xz - 2.*wy,      2.*yz + 2.*wx,      1. - 2.*xx - 2.*yy;
    return out;
  }
  
  SO3 inverse() const
  {
    SO3 q;
    q.arr_(0) = arr_(0);
    q.arr_(1) = -arr_(1);
    q.arr_(2) = -arr_(2);
    q.arr_(3) = -arr_(3);
    return q;
  }
  
  SO3& invert()
  {
    arr_.template block<3,1>(1,0) *= (T)-1.0;
  }
  
  T roll() const
  {
    return atan2(T(2.0)*(w()*x() + y()*z()), T(1.0) - T(2.0)*(x()*x() + y()*y()));
  }

  T pitch() const
  {
    const T val = T(2.0) * (w()*y() - x()*z());

    // hold at 90 degrees if invalid
    if (fabs(val) > T(1.0))
      return copysign(T(1.0), val) * T(M_PI) / T(2.0);
    else
      return asin(val);
  }

  T yaw() const
  {
    return atan2(T(2.0)*(w()*z() + x()*y()), T(1.0) - T(2.0)*(y()*y() + z()*z()));
  }

  Vec3T toEuler() const
  {
    Vec3T out;
    out << roll(), pitch(), yaw();
    return out;
  }
  
  Mat4T qMatLeft()
  {
    Mat4T qL;
    qL << w(), -x(), -y(), -z(),
          x(),  w(), -z(),  y(),
          y(),  z(),  w(), -x(),
          z(), -y(),  x(),  w();
    return qL;
  }
  
  SO3& operator= (const SO3 &q) 
  { 
    arr_ = q.array(); 
  }
  
  template<typename T2>
  SO3 operator* (const SO3<T2> &q) const
  {
    return SO3.fromQuat(qMatLeft() * q.array());
  }
  
  template<typename Tout=T, typename T2>
  Matrix<Tout,3,1> operator* (const Matrix<T2,3,1> &v) const
  {
    Vec3T qv = arr_.template block<3,1>(1,0);
    Matrix<Tout,3,1> t = (Tout)2.0 * v.cross(qv);
    return v - w() * t + t.cross(qv);
  }
  
  template<typename T2>
  SO3& operator*= (const SO3<T2> &q)
  {
    arr_ = qMatLeft() * q.array();
  }
  
  SO3 operator+ (const Vec3T &v) const
  {
    return SO3.fromQuat(qMatLeft() * SO3.Exp(v).array());
  }
  
  SO3& operator+= (const Vec3T &v)
  {
    arr_ = qMatLeft() * SO3.Exp(v).array();
  }
  
  template<typename T2>
  Vec3T operator- (const SO3<T2>& q) const
  {
    return SO3.Log(SO3.fromQuat(q.inverse().qMatLeft() * array()));
  }
  
  static Mat3T hat(const Vec3T &omega)
  {
    Mat3T Omega;
    Omega <<   (T)0, -omega.z(),  omega.y(),
              omega.z(),   (T)0, -omega.x(),
             -omega.y(),  omega.x(),   (T)0;
    return Omega;
  }
  
  static Vec3T vee(const Mat3T &Omega)
  {
    Vec3T omega;
    omega << Omega(2,1), Omega(0,2), Omega(1,0);
    return omega;
  }
  
  static Mat3T log(const SO3 &q)
  {
    return SO3.hat(SO3.Log(q));
  }
  
  static Vec3T Log(const SO3 &q)
  {
    Vec3T qv = q.arr_.template block<3,1>(1,0);
    T     qw = q.w();
    
    T n = qv.norm();
    if (n > 0)
      return 2.0 * qv * atan2(n, qw) / n;
    else
      return Vec3T.Zero();
  }
  
  static SO3 exp(const Mat3T &Omega)
  {
    return SO3.Exp(SO3.vee(Omega));
  }
  
  static SO3 Exp(const Vec3T &omega)
  {
    T th = omega.norm();
    
    SO3 q;
    if (th > 1e-4)
    {
      Vec3T u = omega / th;
      q.arr_(0) = cos(th/2);
      q.arr_.template block<3,1>(1,0) = sin(th/2) * u;
    }
    else
    {
      q.arr_(0) = (T)1.0;
      q.arr_.template block<3,1>(1,0) = omega / 2.0;
      q.normalize();
    }
    
    return q;
  }

};

template<typename T>
inline std::ostream& operator<< (std::ostream& os, const SO3<T>& q)
{
  os << "SO(3): [ " << q.w() << ", " << q.x() << "i, " << q.y() << "j, " << q.z() << "k]";
  return os;
}

typedef SO3<double> SO3d;
