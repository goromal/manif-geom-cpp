#include "SO2.h"
#include "SO3.h"
#include "SE2.h"
#include "SE3.h"

/**
 * @mainpage manif-geom-cpp Library Documentation
 *
 * @section intro_sec Introduction
 * Templated, header-only library implementing the most common manifolds in robotics in 2D and 3D. Operationally very
 * similar to variations on Eigen's `Quaternion<T>` class, but with added
 * chart maps and rules for addition and subtraction on tangent spaces. Meant to be used with nonlinear least-squares
 * solvers like Ceres Solver (http://ceres-solver.org/) which take advantage of templating to implement
 * auto-differentiation on arbitrary mathematical formulations in code.
 *
 * The \f$SO(3)\f$ math is based on [my notes on 3D rotation
 * representations](https://andrewtorgesen.com/notes/Autonomy/Math_Fundamentals/3D_Geometry/Rotations_Robotics_Field_Guide.html).
 *
 * For more information, see
 * [https://andrewtorgesen.com/anixpkgs/cpp/manif-geom-cpp.html](https://andrewtorgesen.com/anixpkgs/cpp/manif-geom-cpp.html).
 *
 * @section install Installation
 * This code is meant to be built as a static library with CMake. It should be compatible with the latest versions of
 * Eigen and Boost (unit test framework only).
 *
 * Install with
 *
 * ```bash
 * mkdir build
 * cd build
 * cmake ..
 * make # or make install
 * ```
 *
 * By default, building will also build and run the unit tests, but this can be turned off with the CMake option
 * `BUILD_TESTS`.
 */
