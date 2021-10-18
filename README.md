# Manifold Geometry C++ Library

Templated, header-only library implementing the [SO(3)](include/SO3.h) and [SE(3)](include/SE3.h) manifolds. Operationally very similar to variations on Eigen's `Quaternion<T>` class, but with added chart maps and rules for addition and subtraction on tangent spaces. Meant to be used with nonlinear least-squares solvers like [Ceres Solver](http://ceres-solver.org/).

