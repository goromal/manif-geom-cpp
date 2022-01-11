#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include "SO3.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestSO3)

BOOST_AUTO_TEST_CASE(TestSO3Action)
{
    static const unsigned int numTests = 1000;
    double quatActNS = 0;
    double rmatActNS = 0;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO3d q = SO3d::random();
        Vector3d v;
        v.setRandom();
        Quaterniond e = Quaterniond::UnitRandom(); // ----

        auto start1 = std::chrono::high_resolution_clock::now();
        Vector3d qv1 = q * v;
        auto stop1 = std::chrono::high_resolution_clock::now();
        quatActNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1-start1).count();

        auto start2 = std::chrono::high_resolution_clock::now();
        Vector3d ev2 = e * v;
        // Vector3d qv2 = q.R() * v;
        auto stop2 = std::chrono::high_resolution_clock::now();
        Vector3d qv2 = q.R() * v; // ----
        rmatActNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop2-start2).count();

        BOOST_CHECK_CLOSE(qv1.x(), qv2.x(), 1e-8);
        BOOST_CHECK_CLOSE(qv1.y(), qv2.y(), 1e-8);
        BOOST_CHECK_CLOSE(qv1.z(), qv2.z(), 1e-8);
    }
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Quaternion action total time: " << quatActNS/1000 << " us" << std::endl;
    std::cout << "Matrix action total time:     " << rmatActNS/1000 << " us" << std::endl;
}

BOOST_AUTO_TEST_CASE(TestInversionAndComposition)
{
    static const unsigned int numTests = 1000;
    double quatInvNS = 0;
    double rmatInvNS = 0;
    double quatCmpNS = 0;
    double rmatCmpNS = 0;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO3d q1 = SO3d::random();
        SO3d q2 = SO3d::random();
        Matrix3d R1 = q1.R();
        Matrix3d R2 = q2.R();
        Quaterniond e1 = Quaterniond::UnitRandom(); // ----
        Quaterniond e2 = Quaterniond::UnitRandom(); // ----

        auto start1 = std::chrono::high_resolution_clock::now();
        SO3d q2inv = q2.inverse();
        auto stop1 = std::chrono::high_resolution_clock::now();
        quatInvNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1-start1).count();

        auto start2 = std::chrono::high_resolution_clock::now();
        Matrix3d R2inv = R2.transpose();
        auto stop2 = std::chrono::high_resolution_clock::now();
        rmatInvNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop2-start2).count();

        auto start3 = std::chrono::high_resolution_clock::now();
        SO3d q1p = q1 * q2 * q2inv;
        auto stop3 = std::chrono::high_resolution_clock::now();
        quatCmpNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop3-start3).count();

        auto start4 = std::chrono::high_resolution_clock::now();
        // Matrix3d R1p = R1 * R2 * R2inv;
        Quaterniond e3 = e1 * e2 * e2;
        auto stop4 = std::chrono::high_resolution_clock::now();
        rmatCmpNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop4-start4).count();
        
        BOOST_CHECK_CLOSE(q1.w(), q1p.w(), 1e-8);
        BOOST_CHECK_CLOSE(q1.x(), q1p.x(), 1e-8);
        BOOST_CHECK_CLOSE(q1.y(), q1p.y(), 1e-8);
        BOOST_CHECK_CLOSE(q1.z(), q1p.z(), 1e-8);
    }
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Quaternion inversion total time: " << quatInvNS/1000 << " us" << std::endl;
    std::cout << "Matrix inversion total time:     " << rmatInvNS/1000 << " us" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Quaternion multiply total time:  " << quatCmpNS/1000 << " us" << std::endl;
    std::cout << "Matrix multiply total time:      " << rmatCmpNS/1000 << " us" << std::endl;
}

BOOST_AUTO_TEST_CASE(TestEulerConvs)
{
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        Vector3d euler;
        euler.setRandom();
        euler *= M_PI;
        SO3d q = SO3d::fromEuler(euler.x(), euler.y(), euler.z());
        SO3d q2 = SO3d::fromEuler(q.roll(), q.pitch(), q.yaw());

        BOOST_CHECK_CLOSE(q.w(), q2.w(), 1e-8);
        BOOST_CHECK_CLOSE(q.x(), q2.x(), 1e-8);
        BOOST_CHECK_CLOSE(q.y(), q2.y(), 1e-8);
        BOOST_CHECK_CLOSE(q.z(), q2.z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO3d q1 = SO3d::random();
        Vector3d q12;
        q12.setRandom();
        SO3d q2 = q1 + q12;
        Vector3d q12p = q2 - q1;
        BOOST_CHECK_CLOSE(q12.x(), q12p.x(), 1e-8);
        BOOST_CHECK_CLOSE(q12.y(), q12p.y(), 1e-8);
        BOOST_CHECK_CLOSE(q12.z(), q12p.z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestConstructors)
{
    Vector4d q_vec(1., 1., 0., 0.);
    SO3d q(q_vec);
    BOOST_CHECK_CLOSE(q.w(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q.x(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q.y(), 0.0, 1e-8);
    BOOST_CHECK_CLOSE(q.z(), 0.0, 1e-8);

    double q_arr[4];
    q_arr[0] = 1.;
    q_arr[1] = 1.;
    q_arr[2] = 0.;
    q_arr[3] = 0.;
    SO3d q2(q_arr);
    BOOST_CHECK_CLOSE(q2.w(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q2.x(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q2.y(), 0.0, 1e-8);
    BOOST_CHECK_CLOSE(q2.z(), 0.0, 1e-8);
}

BOOST_AUTO_TEST_CASE(TestMutableArray)
{
    SO3d q = SO3d::identity();
    Vector4d q_arr = q.array();
    q_arr(0) = 2.;
    BOOST_CHECK_CLOSE(q_arr(0), 2., 1e-8);
    BOOST_CHECK_CLOSE(q.w(), 1.0, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()