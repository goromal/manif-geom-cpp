#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include "SO2.h"

using namespace Eigen;

BOOST_AUTO_TEST_SUITE(TestSO2)

BOOST_AUTO_TEST_CASE(TestAction)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    double so2ActNS = 0;
    double so3ActNS = 0;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO2d so2 = SO2d::random();
        Vector2d v2;
        v2.setRandom();
        Vector3d v3;
        v3 << v2.x(), v2.y(), 0;
        SO3d so3 = SO3d::fromEuler(0, 0, so2.angle());

        auto start1 = std::chrono::high_resolution_clock::now();
        Vector2d qv2 = so2 * v2;
        auto stop1 = std::chrono::high_resolution_clock::now();
        so2ActNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1-start1).count();

        auto start2 = std::chrono::high_resolution_clock::now();
        Vector3d qv3 = so3 * v3;
        auto stop2 = std::chrono::high_resolution_clock::now();
        so3ActNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop2-start2).count();

        BOOST_CHECK_CLOSE(qv2.x(), qv3.x(), 1e-8);
        BOOST_CHECK_CLOSE(qv2.y(), qv3.y(), 1e-8);
    }
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "SO2d action average time: " << so2ActNS/1000 << " ns" << std::endl;
    std::cout << "SO3d action average time: " << so3ActNS/1000 << " ns" << std::endl;
}

BOOST_AUTO_TEST_CASE(TestInversionAndComposition)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    double so2CmpNS = 0;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO2d q1 = SO2d::random();
        SO2d q2 = SO2d::random();

        SO2d q2inv = q2.inverse();

        auto start1 = std::chrono::high_resolution_clock::now();
        SO2d q1p = q1 * q2 * q2inv;
        auto stop1 = std::chrono::high_resolution_clock::now();
        so2CmpNS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1-start1).count();

        BOOST_CHECK_CLOSE(q1.w(), q1p.w(), 1e-8);
        BOOST_CHECK_CLOSE(q1.x(), q1p.x(), 1e-8);
    }
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "SO2d multiply (3x) average time: " << so2CmpNS/1000 << " ns" << std::endl;
}

BOOST_AUTO_TEST_CASE(TestAngleConvs)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        Vector1d angle;
        angle.setRandom();
        angle *= M_PI;
        SO2d q = SO2d::fromAngle(angle.x());
        SO2d q2 = SO2d::fromEuler(q.angle());

        BOOST_CHECK_CLOSE(q.w(), q2.w(), 1e-8);
        BOOST_CHECK_CLOSE(q.x(), q2.x(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO2d q1 = SO2d::random();
        Vector1d q12;
        q12.setRandom();
        SO2d q2 = q1 + q12;
        Vector1d q12p = q2 - q1;
        BOOST_CHECK_CLOSE(q12.x(), q12p.x(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestChartMaps)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO2d q = SO2d::random();
        Vector1d w;
        w.setRandom();
        Vector1d qLog = SO2d::Log(q);
        SO2d q2 = SO2d::Exp(qLog);
        BOOST_CHECK_CLOSE(q.w(), q2.w(), 1e-8);
        BOOST_CHECK_CLOSE(q.x(), q2.x(), 1e-8);
        SO2d wExp = SO2d::Exp(w);
        Vector1d w2 = SO2d::Log(wExp);
        BOOST_CHECK_CLOSE(w.x(), w2.x(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestConstructors)
{
    Vector2d q_vec(1., 0.);
    SO2d q(q_vec);
    BOOST_CHECK_CLOSE(q.w(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q.x(), 0.0, 1e-8);

    double q_arr[2];
    q_arr[0] = 1.;
    q_arr[1] = 0.;
    SO2d q2(q_arr);
    BOOST_CHECK_CLOSE(q2.w(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(q2.x(), 0.0, 1e-8);
}

BOOST_AUTO_TEST_CASE(TestMutableArray)
{
    SO2d q = SO2d::identity();
    Vector2d q_arr = q.array();
    q_arr(0) = 2.;
    BOOST_CHECK_CLOSE(q_arr(0), 2., 1e-8);
    BOOST_CHECK_CLOSE(q.w(), 1.0, 1e-8);
}

BOOST_AUTO_TEST_CASE(TestScaling)
{
    srand(444444);
    SO2d qI = SO2d::identity();
    SO2d qIs = 5.0 * qI;
    BOOST_CHECK_CLOSE(qIs.w(), qI.w(), 1e-8);
    BOOST_CHECK_CLOSE(qIs.x(), qI.x(), 1e-8);
    SO2d qr = SO2d::random();
    SO2d qr2 = qr * 0.2; // if scale is too big, then the rotation will
                         // wrap around the circle, resulting in a reversed
                         // or truncated tangent vector which can't be inverted
                         // through scalar division
    SO2d qr3 = qr2 / 0.2;
    BOOST_CHECK_CLOSE(qr.w(), qr3.w(), 1e-8);
    BOOST_CHECK_CLOSE(qr.x(), qr3.x(), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
