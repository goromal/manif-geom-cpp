# Manifold Geometry C++ Library

![example workflow](https://github.com/goromal/manif-geom-cpp/actions/workflows/test.yml/badge.svg)

Templated, header-only library implementing the [SO(3)](include/SO3.h) and [SE(3)](include/SE3.h) manifolds, along with their 2D corollaries. Operationally very similar to variations on Eigen's `Quaternion<T>` class, but with added chart maps and rules for addition and subtraction on tangent spaces. Meant to be used with nonlinear least-squares solvers like [Ceres Solver](http://ceres-solver.org/) which take advantage of templating to implement auto-differentiation on arbitrary mathematical formulations in code.

The SO3 math is based on [my notes on 3D rotation representations](https://notes.andrewtorgesen.com/doku.php?id=public:autonomy:math:3d-geometry:implementing-rotations).

## Building / Installing

This library is built with CMake. Most recently tested with the following dependency versions:

- Eigen 3.4.0
- Boost 1.79.0 (for unit test framework)

```bash
mkdir build
cd build
cmake ..
make # or make install
```

By default, building will build and run the unit tests, but this can be turned off with the CMake option `BUILD_TESTS`.
