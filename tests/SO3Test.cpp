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
    double quatActMS = 0;
    double rmatActMS = 0;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SO3d q = SO3d::random();
        Vector3d v;
        v.setRandom();

        auto start1 = std::chrono::high_resolution_clock::now();
        Vector3d qv1 = q * v;
        auto stop1 = std::chrono::high_resolution_clock::now();
        quatActMS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1-start1).count();

        auto start2 = std::chrono::high_resolution_clock::now();
        Vector3d qv2 = q.R() * v;
        auto stop2 = std::chrono::high_resolution_clock::now();
        rmatActMS += std::chrono::duration_cast<std::chrono::nanoseconds>(stop2-start2).count();

        BOOST_CHECK_CLOSE(qv1.x(), qv2.x(), 1e-8);
        BOOST_CHECK_CLOSE(qv1.y(), qv2.y(), 1e-8);
        BOOST_CHECK_CLOSE(qv1.z(), qv2.z(), 1e-8);
    }
    std::cout << "Quaternion action total time: " << quatActMS/1000 << " us" << std::endl;
    std::cout << "Matrix action total time:     " << rmatActMS/1000 << " us" << std::endl;
}

BOOST_AUTO_TEST_CASE(TestEulerConvs)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestRotMatConvs)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestUVInit)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestInversion)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestComposition)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestChartMaps)
{
    // TODO   
}

BOOST_AUTO_TEST_SUITE_END()