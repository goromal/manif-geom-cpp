#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "SO3.h"
#include <math.h>
#include <iostream>

using namespace Eigen;

template<typename T>
class SE3
{

private:
  typedef Matrix<T,3,1> Vec3T;
  typedef Matrix<T,4,1> Vec4T;
  typedef Matrix<T,6,1> Vec6T;
  typedef Matrix<T,7,1> Vec7T;
  typedef Matrix<T,3,3> Mat3T;
  typedef Matrix<T,4,4> Mat4T;
  typedef Matrix<T,6,6> Mat6T;
  
  T buf_[7];

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Map<Vec7T> arr_;
  Map<Vec3T> t_;
  SO3<T> q_;
  
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
  
  static SE3 fromH(const Mat4T &m)
  {
    SE3 x;
    x.t_ = m.template block<3,1>(0,3);
    x.q_ = SO3<T>::fromR(m.template block<3,3>(0,0));
    return x;
  }
  
  static SE3 fromVecAndQuat(const T tx, const T ty, const T tz, const T qw, const T qx, const T qy, const T qz)
  {
    SE3 x;
    x.arr_ << tx, ty, tz, qw, qx, qy, qz;
    return x;
  }
  
  static SE3 fromVecAndQuat(const Vec3T &tvec, const Vec4T &qvec)
  {
    SE3 x;
    x.arr_ << tvec(0), tvec(1), tvec(2), qvec(0), qvec(1), qvec(2), qvec(3);
    return x;
  }
  
  static SE3 fromVecAndQuat(const Vec3T &tvec, const Quaternion<T> quat)
  {
    SE3 x;
    x.arr_ << tvec(0), tvec(1), tvec(2), quat.w(), quat.x(), quat.y(), quat.z();
    return x;
  }
  
  SE3() :
    arr_(buf_),
    t_(arr_.data()),
    q_(arr_.data()+3)
  {}

  SE3(const Ref<const Vec7T> arr) :
    arr_(const_cast<T*>(arr.data())),
    t_(arr_.data()),
    q_(arr_.data()+3)
  {}

  SE3(const SE3& x) : 
    arr_(buf_),
    t_(arr_.data()),
    q_(arr_.data()+3)
  {
    arr_ = x.arr_;
  }

  SE3(const T* data) :
    arr_(const_cast<T*>(data)),
    t_(arr_.data()),
    q_(arr_.data()+3)
  {}

  inline T* data() { return arr_.data(); }
  inline const T* data() const { return arr_.data();}
  inline T& operator[] (int i) {return arr_[i];}
  inline Vec3T& t() const { return t_; }
  inline SO3<T>& q() const { return q_; }
  inline const Vec7T array() const { return arr_;}
  
  Mat4T H() const
  {
    Mat4T out;
    out << q().R(), t(),
           (T)0, (T)0, (T)0, (T)1;
    return out;
  }
  
  SE3 inverse() const
  {
    SE3 x;
    SO3<T> q_inv = q().inverse();
    x.t() = -q_inv * t();
    x.q() = q_inv;
    return x;
  }
  
  SE3& invert()
  {
    q().invert();
    t() = -q() * t();
  }
  
  SE3& operator= (const SE3 &x) 
  { 
    arr_ = x.array(); 
  }
  
  template<typename T2>
  SE3 operator* (const SE3<T2> &x) const
  {
    SE3 x_out;
    x_out.t() = t() + q() * x.t();
    x_out.q() = q() * x.q();
    return x_out;
  }
  
  template<typename Tout=T, typename T2>
  Matrix<Tout,3,1> operator* (const Matrix<T2,3,1> &v) const
  {
    return q() * v + t();
  }
  
  template<typename T2>
  SE3& operator*= (const SE3<T2> &x)
  {
    t() = t() + q() * x.t();
    q() = q() * x.q();
  }
  
  SE3 operator+ (const Vec6T &v) const
  {
    return *this * SE3::Exp(v);
  }
  
  SE3& operator+= (const Vec6T &v)
  {
    SE3 x = SE3::Exp(v);
    t() = t() + q() * x.t();
    q() = q() * x.q();
  }
  
  template<typename T2>
  Vec6T operator- (const SE3<T2>& x) const
  {
    return SE3::Log(x.inverse() * *this);
  }
  
  static Mat6T hat(const Vec6T &omega)
  {
    Mat6T Omega;
    Omega << SO3<T>::hat(omega.template block<3,1>(3,0)), omega.template block<3,1>(0,0),
             (T)0, (T)0, (T)0, (T)0;
    return Omega;
  }
  
  static Vec6T vee(const Mat6T &Omega)
  {
    Vec6T omega;
    omega << Omega.template block<3,1>(0,3), SO3<T>::vee(Omega.template block<3,3>(0,0));
    return omega;
  }
  
  static Mat6T log(const SE3 &x)
  {
    return SE3::hat(SE3::Log(x));
  }
  
  static Vec6T Log(const SE3 &x)
  {
    Vec3T w = SO3<T>::Log(x.q());
    T    th = w.norm();
    Mat3T W = SO3<T>::hat(w);
    
    Mat3T leftJacobianInverse;
    if (th > 0)
    {
      T a = sin(th)/th;
      T b = ((T)1. - cos(th))/(th*th);
      T c = ((T)1. - a)/(th*th);
      T e = (b - 2.*c)/(2.*a);
      leftJacobianInverse = Mat3T::Identity() - 0.5*W + e*(W*W);
    }
    else
    {
      leftJacobianInverse = Mat3T::Identity();
    }
    
    Vec6T out;
    out << leftJacobianInverse * x.t(), w;
    
    return out;
  }
  
  static SE3 exp(const Mat6T &Omega)
  {
    return SE3::Exp(SE3::vee(Omega));
  }
  
  static SE3 Exp(const Vec6T &omega)
  {
    Vec3T rho = omega.template block<3,1>(0,0);
    Vec3T w = omega.template block<3,1>(3,0);
    SO3<T> q = SO3<T>::Exp(w);
    Mat3T W = SO3<T>::hat(w);
    T th = w.norm();
  
    Mat3T leftJacobian;
    if (th > 0)
    {
      T a = sin(th)/th;
      T b = ((T)1. - cos(th))/(th*th);
      T c = ((T)1. - a)/(th*th);
      leftJacobian = a * Mat3T::Identity() + b*W + c*(w*w.transpose());
    }
    else
    {
      leftJacobian = Mat3T::Identity();
    }
    
    SE3 x;
    x.t() = leftJacobian * rho;
    x.q() = q;
    
    return x;
  }

};

template<typename T>
inline std::ostream& operator<< (std::ostream& os, const SE3<T>& x)
{
  os << "SE(3): [ " << x.t().x() << "i, " << x.t().y() << "j, " << x.t().z() << "k] [ " << x.q().w() << ", " << x.q().x() << "i, " << x.q().y() << "j, " << x.q().z() << "k]";
  return os;
}

typedef SE3<double> SE3d;
