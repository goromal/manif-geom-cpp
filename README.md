# Manifold Geometry C++ Library

![example workflow](https://github.com/goromal/manif-geom-cpp/actions/workflows/test.yml/badge.svg)

Templated, header-only library implementing the [SO(3)](include/SO3.h) and [SE(3)](include/SE3.h) manifolds, along with their 2D corollaries. Operationally very similar to variations on Eigen's `Quaternion<T>` class, but with added chart maps and rules for addition and subtraction on tangent spaces. Meant to be used with nonlinear least-squares solvers like [Ceres Solver](http://ceres-solver.org/) which take advantage of templating to implement auto-differentiation on arbitrary mathematical formulations in code.

The SO3 math is based on [my notes on 3D rotation representations](https://notes.andrewtorgesen.com/doku.php?id=public:implementing-rotations).
