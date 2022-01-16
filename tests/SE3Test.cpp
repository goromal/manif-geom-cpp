#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <chrono>
#include "SE3.h"

using namespace Eigen;
typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,7,1> Vector7d;

BOOST_AUTO_TEST_SUITE(TestSE3)

BOOST_AUTO_TEST_CASE(TestSE3Action)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(TestInversionAndComposition)
{
    // TODO   
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE3d x1 = SE3d::random();
        Vector6d x12;
        x12.setRandom();
        SE3d x2 = x1 + x12;
        Vector6d x12p = x2 - x1;
        BOOST_CHECK_CLOSE(x12(0), x12p(0), 1e-8);
        BOOST_CHECK_CLOSE(x12(1), x12p(1), 1e-8);
        BOOST_CHECK_CLOSE(x12(2), x12p(2), 1e-8);
        BOOST_CHECK_CLOSE(x12(3), x12p(3), 1e-8);
        BOOST_CHECK_CLOSE(x12(4), x12p(4), 1e-8);
        BOOST_CHECK_CLOSE(x12(5), x12p(5), 1e-8);
    }  
}

BOOST_AUTO_TEST_CASE(TestConstructors)
{
    Vector7d x_vec(1., 1., 0., 1., 0., 0., 0.);
    SE3d x(x_vec);
    BOOST_CHECK_CLOSE(x.t().x(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x.t().y(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x.t().z(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x.q().w(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x.q().x(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x.q().y(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x.q().z(), 0., 1e-8);

    double x_arr[7];
    x_arr[0] = 1.;
    x_arr[1] = 1.;
    x_arr[2] = 0.;
    x_arr[3] = 1.;
    x_arr[4] = 0.;
    x_arr[5] = 0.;
    x_arr[6] = 0.;
    SE3d x2(x_arr);
    BOOST_CHECK_CLOSE(x2.t().x(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x2.t().y(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x2.t().z(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().w(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().x(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().y(), 0., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().z(), 0., 1e-8);
}

BOOST_AUTO_TEST_CASE(TestMutableArray)
{
    SE3d x = SE3d::identity();
    Vector7d x_arr = x.array();
    x_arr(0) = 2.;
    BOOST_CHECK_CLOSE(x_arr(0), 2., 1e-8);
    BOOST_CHECK_CLOSE(x.t().x(), 0., 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()