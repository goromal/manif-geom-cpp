#include "SE2.h"
#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <chrono>

using namespace Eigen;
typedef Matrix<double, 1, 1> Vector1d;
typedef Matrix<double, 2, 1> Vector2d;
typedef Matrix<double, 4, 1> Vector4d;

BOOST_AUTO_TEST_SUITE(TestSE2)

BOOST_AUTO_TEST_CASE(TestAction)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE2d     x    = SE2d::random();
        SE2d     xinv = x.inverse();
        Vector2d v;
        v.setRandom();
        Vector2d v2 = xinv * (x * v);
        BOOST_CHECK_CLOSE(v.x(), v2.x(), 1e-8);
        BOOST_CHECK_CLOSE(v.y(), v2.y(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestInversionAndComposition)
{
    srand(444444);
    static const unsigned int numTests = 1000;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE2d x1    = SE2d::random();
        SE2d x2    = SE2d::random();
        SE2d x2inv = x2.inverse();
        SE2d x1p   = x1 * x2 * x2inv;
        BOOST_CHECK_CLOSE(x1.t().x(), x1p.t().x(), 1e-8);
        BOOST_CHECK_CLOSE(x1.t().y(), x1p.t().y(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().w(), x1p.q().w(), 1e-8);
        BOOST_CHECK_CLOSE(x1.q().x(), x1p.q().x(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestPlusMinus)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE2d     x1 = SE2d::random();
        Vector3d x12;
        x12.setRandom();
        SE2d     x2   = x1 + x12;
        Vector3d x12p = x2 - x1;
        BOOST_CHECK_CLOSE(x12(0), x12p(0), 1e-8);
        BOOST_CHECK_CLOSE(x12(1), x12p(1), 1e-8);
        BOOST_CHECK_CLOSE(x12(2), x12p(2), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestChartMaps)
{
    srand(444444);
    static const unsigned int numTests = 50;
    for (unsigned int i = 0; i < numTests; i++)
    {
        SE2d     x = SE2d::random();
        Vector3d w;
        w.setRandom();
        Vector3d xLog = SE2d::Log(x);
        SE2d     x2   = SE2d::Exp(xLog);
        BOOST_CHECK_CLOSE(x.t().x(), x2.t().x(), 1e-8);
        BOOST_CHECK_CLOSE(x.t().y(), x2.t().y(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().w(), x2.q().w(), 1e-8);
        BOOST_CHECK_CLOSE(x.q().x(), x2.q().x(), 1e-8);
        SE2d     wExp = SE2d::Exp(w);
        Vector3d w2   = SE2d::Log(wExp);
        BOOST_CHECK_CLOSE(w.x(), w2.x(), 1e-8);
        BOOST_CHECK_CLOSE(w.y(), w2.y(), 1e-8);
        BOOST_CHECK_CLOSE(w.z(), w2.z(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(TestConstructors)
{
    Vector4d x_vec(1., -2., 1., 0.);
    SE2d     x(x_vec);
    BOOST_CHECK_CLOSE(x.t().x(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x.t().y(), -2., 1e-8);
    BOOST_CHECK_CLOSE(x.q().w(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x.q().x(), 0., 1e-8);

    double x_arr[4];
    x_arr[0] = 1.;
    x_arr[1] = -2.;
    x_arr[2] = 1.;
    x_arr[3] = 0.;
    SE2d x2(x_arr);
    BOOST_CHECK_CLOSE(x2.t().x(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x2.t().y(), -2., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().w(), 1., 1e-8);
    BOOST_CHECK_CLOSE(x2.q().x(), 0., 1e-8);
}

BOOST_AUTO_TEST_CASE(TestMutableArray)
{
    SE2d     x     = SE2d::identity();
    Vector4d x_arr = x.array();
    x_arr(0)       = 2.;
    BOOST_CHECK_CLOSE(x_arr(0), 2., 1e-8);
    BOOST_CHECK_CLOSE(x.t().x(), 0., 1e-8);
}

BOOST_AUTO_TEST_CASE(TestScaling)
{
    srand(444444);
    SE2d xI  = SE2d::identity();
    SE2d xIs = 5.0 * xI;
    BOOST_CHECK_CLOSE(xIs.t().x(), xI.t().x(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.t().y(), xI.t().y(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().w(), xI.q().w(), 1e-8);
    BOOST_CHECK_CLOSE(xIs.q().x(), xI.q().x(), 1e-8);
    SE2d xr  = SE2d::random();
    SE2d xr2 = xr * 0.2;
    SE2d xr3 = xr2 / 0.2;
    BOOST_CHECK_CLOSE(xr.t().x(), xr3.t().x(), 1e-8);
    BOOST_CHECK_CLOSE(xr.t().y(), xr3.t().y(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().w(), xr3.q().w(), 1e-8);
    BOOST_CHECK_CLOSE(xr.q().x(), xr3.q().x(), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
