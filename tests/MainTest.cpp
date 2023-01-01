#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ManifGeomCppTests

#include <boost/test/unit_test.hpp>
#include <iostream>

/*! Global testing definitions. */
struct GlobalFixture
{
    GlobalFixture()
    {
        // common function calls here
    }

    ~GlobalFixture() {}
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);
