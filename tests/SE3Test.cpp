#include "SE3.h"
#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <chrono>

using namespace Eigen;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 7, 1> Vector7d;

BOOST_AUTO_TEST_SUITE(TestSE3)

BOOST_AUTO_TEST_CASE(TestAction)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE3d     x    = SE3d::random();
        SE3d     xinv = x.inverse();
        Vector3d v;
        v.setRandom();
        Vector3d v2 = xinv * (x * v);
        BOOST_CHECK_CLOSE(v.x(), v2.x(), 1e-8);
        BOOST_CHECK_CLOSE(v.y(), v2.y(), 1e-8);
        BOOST_CHECK_CLOSE(v.z(), v2.z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestInversionAndComposition)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE3d x1    = SE3d::random();
        SE3d x2    = SE3d::random();
        SE3d x2inv = x2.inverse();
        SE3d x1p   = x1 * x2 * x2inv;
        BOOST_CHECK_CLOSE(x1.t().x(), x1p.t().x(), 1e-8);
        BOOST_CHECK_CLOSE(x1.t().y(), x1p.t().y(), 1e-8);
        BOOST_CHECK_CLOSE(x1.t().z(), x1p.t().z(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().w(), x1p.q().w(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().x(), x1p.q().x(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().y(), x1p.q().y(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().z(), x1p.q().z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE3d     x1 = SE3d::random();
        Vector6d x12;
        x12.setRandom();
        SE3d     x2   = x1 + x12;
        Vector6d x12p = x2 - x1;
        BOOST_CHECK_CLOSE(x12(0), x12p(0), 1e-8);
        BOOST_CHECK_CLOSE(x12(1), x12p(1), 1e-8);
        BOOST_CHECK_CLOSE(x12(2), x12p(2), 1e-8);
        BOOST_CHECK_CLOSE(x12(3), x12p(3), 1e-8);
        BOOST_CHECK_CLOSE(x12(4), x12p(4), 1e-8);
        BOOST_CHECK_CLOSE(x12(5), x12p(5), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestChartMaps)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE3d                 x = SE3d::random();
        Matrix<double, 6, 1> w;
        w.setRandom();
        Matrix<double, 6, 1> xLog = SE3d::Log(x);
        SE3d                 x2   = SE3d::Exp(xLog);
        BOOST_CHECK_CLOSE(x.t().x(), x2.t().x(), 1e-8);
        BOOST_CHECK_CLOSE(x.t().y(), x2.t().y(), 1e-8);
        BOOST_CHECK_CLOSE(x.t().z(), x2.t().z(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().w(), x2.q().w(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().x(), x2.q().x(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().y(), x2.q().y(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().z(), x2.q().z(), 1e-8);
        SE3d                 wExp = SE3d::Exp(w);
        Matrix<double, 6, 1> w2   = SE3d::Log(wExp);
        BOOST_CHECK_CLOSE(w.x(), w2.x(), 1e-8);
        BOOST_CHECK_CLOSE(w.y(), w2.y(), 1e-8);
        BOOST_CHECK_CLOSE(w.z(), w2.z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestConstructors)
{
    Vector7d x_vec(1., 1., 0., 1., 0., 0., 0.);
    SE3d     x(x_vec);
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
    SE3d     x     = SE3d::identity();
    Vector7d x_arr = x.array();
    x_arr(0)       = 2.;
    BOOST_CHECK_CLOSE(x_arr(0), 2., 1e-8);
    BOOST_CHECK_CLOSE(x.t().x(), 0., 1e-8);
}

BOOST_AUTO_TEST_CASE(TestScaling)
{
    srand(444444);
    SE3d xI  = SE3d::identity();
    SE3d xIs = 5.0 * xI;
    BOOST_CHECK_CLOSE(xIs.t().x(), xI.t().x(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.t().y(), xI.t().y(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.t().z(), xI.t().z(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().w(), xI.q().w(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().x(), xI.q().x(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().y(), xI.q().y(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().z(), xI.q().z(), 1e-8);
    SE3d xr  = SE3d::random();
    SE3d xr2 = xr * 0.2;
    SE3d xr3 = xr2 / 0.2;
    BOOST_CHECK_CLOSE(xr.t().x(), xr3.t().x(), 1e-8);
    BOOST_CHECK_CLOSE(xr.t().y(), xr3.t().y(), 1e-8);
    BOOST_CHECK_CLOSE(xr.t().z(), xr3.t().z(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().w(), xr3.q().w(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().x(), xr3.q().x(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().y(), xr3.q().y(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().z(), xr3.q().z(), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
