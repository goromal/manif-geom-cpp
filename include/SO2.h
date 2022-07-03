#pragma once

#include <Eigen/Core>
#include <math.h>
#include <iostream>

using namespace Eigen;

template<typename T>
class SO2
{
private:
    typedef Matrix<T,1,1> Vec1T;
    typedef Matrix<T,2,1> Vec2T;
    typedef Matrix<T,2,2> Mat2T;
    
    T buf_[2];
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Map<Vec2T> arr_;
    
    static SO2 random()
    {
        SO2 q;
        q.arr_.setRandom();
        q.normalize();
        return q;
    } 
    
    static SO2 identity()
    {
        SO2 q;
        q.arr_ << 1., 0.;
        return q;
    }
    
    static SO2 fromAngle(const T angle)
    {
        SO2 q;
        q.arr_ << cos(angle), sin(angle);
        return q;
    }
    
    static SO2 fromR(const Mat2T &m)
    {
        SO2 q;
        q.arr_ << m(0,0), m(1,0);
        return q;
    }
    
    static SO2 fromTwoUnitVectors(const Vec2T &u, const Vec2T &v)
    {
        SO2 q;
        T d = u.dot(v);
        
        if (d < (T)0.99999999 && d > (T)-0.99999999)
        {
            q.arr_ << d, sqrt(1. - d * d);
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
    
    static SO2 fromComplex(const T qw, const T qx)
    {
        SO2 q;
        q.arr_ << qw, qx;
        return q;
    }
    
    static SO2 fromComplex(const Vec2T &qvec)
    {
        SO2 q;
        q.arr_ << qvec(0), qvec(1);
        return q;
    }
    
    SO2()
      : arr_(buf_)
    {}
    
    SO2(const Ref<const Vec2T>& arr)
      : arr_(buf_)
    {
        arr_ = arr;
    }
    
    SO2(const SO2& q)
      : arr_(buf_)
    {
        arr_ = q.arr_;
    }
    
    SO2(const T* data)
      : arr_(const_cast<T*>(data))
    {}
    
    inline T& operator[] (int i) {return arr_[i];}
    inline const T& w() const { return arr_(0); }
    inline const T& x() const { return arr_(1); }
    inline T& w() { return arr_(0); }
    inline T& x() { return arr_(1); }
    inline const Vec2T elements() const { return arr_; }
    inline Vec2T array() const { return arr_; }
    inline T* data() { return arr_.data(); }
    inline const T* data() const { return arr_.data(); }
    
    SO2 copy() const
    {
        SO2 tmp;
        tmp.arr_ = arr_;
        return tmp;   
    }
    
    void normalize()
    {
        arr_ /= arr_.norm();
    }
    
    SO2 normalized()
    {
        SO2 tmp = copy();
        tmp.normalize();
        return tmp;
    }
    
    Mat2T R() const
    {
        T costh = w();
        T sinth = x();
        Mat2T out;
        out << costh, -sinth,
               sinth, costh;
        return out; 
    }
    
    SO2 inverse() const
    {
        SO2 q;
        q.arr_(0) = arr_(0);
        q.arr_(1) = -arr_(1);
        return q;
    }
    
    SO2& invert()
    {
        arr_(1) *= (T)-1.0;
        return *this;
    }
    
    T angle() const
    {
        return atan2(x(), w());
    }
    
    template <typename Tout=T, typename T2>
    SO2<Tout> otimes(const SO2<T2>& q) const
    {
        SO2<Tout> qout;
        qout.arr_ <<  w() * q.w() - x() * q.x(),
                      w() * q.x() + x() * q.w();
        return qout;
    }

    template<typename Tout=T, typename T2>
    SO2<Tout> oplus(const Matrix<T2,1,1> &delta) const
    {
        return otimes<Tout, T2>(SO2<T2>::Exp(delta));
    }

    template<typename Tout=T, typename T2>
    Matrix<Tout,1,1> ominus(const SO2<T2> &q) const
    {
        SO2<Tout> dq = q.inverse().template otimes<Tout>(*this);
        return SO2<Tout>::Log(dq);
    }

    SO2& operator= (const SO2 &q) 
    { 
        arr_ = q.elements(); 
        return *this;
    }
  
    template<typename T2>
    SO2 operator* (const SO2<T2> &q) const
    {
        return otimes(q);
    }

    template<typename T2>
    SO2& operator*= (const SO2<T2> &q)
    {
        arr_ = otimes(q).elements();
        return *this;
    }

    SO2& operator*= (const double &s)
    {
        arr_ = SO2::Exp(s * SO2::Log(*this)).elements();
        return *this;
    }

    SO2& operator/= (const double &s)
    {
        arr_ = SO2::Exp(SO2::Log(*this) / s).elements();
        return *this;
    }

    SO2 operator/ (const double &s) const
    {
        SO2 qs;
        qs.arr_ = SO2::Exp(SO2::Log(*this) / s).elements();
        return qs;
    }
  
    template<typename Tout=T, typename T2>
    Matrix<Tout,2,1> operator* (const Matrix<T2,2,1> &v) const
    {
        Vec2T out;
        out << w() * v.x() - x() * v.y(),
               w() * v.y() + x() * v.x();
        return out;
    }
  
    Vec2T operator* (const Vec2T &v) const
    {
        Vec2T out;
        out << w() * v.x() - x() * v.y(),
               w() * v.y() + x() * v.x();
        return out;
    }
  
    SO2 operator+ (const Vec1T &v) const
    {
        return oplus(v);
    }
  
    SO2& operator+= (const Vec1T &v)
    {
        arr_ = oplus(v).elements();
        return *this;
    }
  
    template<typename T2>
    Vec1T operator- (const SO2<T2>& q) const
    {
        return ominus(q);
    }
  
    static Mat2T hat(const Vec1T &omega)
    {
        Mat2T Omega;
        Omega <<   (T)0., -omega.x(),
               omega.x(),      (T)0.;
        return Omega;
    }

    static Vec1T vee(const Mat2T &Omega)
    {
        Vec1T omega;
        omega << Omega(1,0);
        return omega;
    }
  
    static Mat2T log(const SO2 &q)
    {
        return hat(SO2::Log(q));
    }

    static Vec1T Log(const SO2 &q)
    {
        Vec1T out;
        out << q.angle();
        return out;
    }
  
    static SO2 exp(const Mat2T &Omega)
    {
        return SO2::Exp(vee(Omega));
    }

    static SO2 Exp(const Vec1T &omega)
    {
        return SO2::fromAngle(omega.x());
    }
  
    template<typename T2>
    SO2<T2> cast() const
    {
        SO2<T2> q;
        q.arr_ = arr_.template cast<T2>();
        return q;
    }
};

template<typename T>
SO2<T> operator* (const double &l, const SO2<T> &r)
{
    SO2<T> lr;
    lr.arr_ = SO2<T>::Exp(l * SO2<T>::Log(r)).elements();
    return lr;
}

template<typename T>
SO2<T> operator* (const SO2<T> &l, const double &r)
{
    SO2<T> lr;
    lr.arr_ = SO2<T>::Exp(r * SO2<T>::Log(l)).elements();
    return lr;
}

template<typename T>
inline std::ostream& operator<< (std::ostream& os, const SO2<T> &q)
{
    os << "SO(2): [ " << q.w() << ", " << q.x() << "i ]";
    return os;
}

typedef SO2<double> SO2d;
