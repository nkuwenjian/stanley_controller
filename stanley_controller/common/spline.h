// =============================================================================
// Copyright (C) 2022 Xiaomi EV Company Limited. All rights reserved.
//
// XIAOMI CONFIDENTIAL
//
// Xiaomi EV Company Limited retains all intellectual property and proprietary
// rights in and to this software and related documentation and any
// modifications thereto ("Software"). Any reproduction, disclosure or
// distribution of this Software is strictly prohibited.
//
// Permitted Users may only access, use or edit the Software based on their work
// duties. Permitted Users shall refer to the current employees in Xiaomi EV
// Autonomous Driving Group who have entered into Non-disclosure Agreement with
// Xiaomi, and have been duly authorized to access the Software.
// =============================================================================

#pragma once

#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "glog/logging.h"

namespace stanley_controller {

/*
 * Solves a tridiagonal system of equations.
 *
 * The arguments lower, diag, upper and rhs define the equations and must
 * be of the same length (say n). The i-th element of each vector specifies
 * the i-th linear equation, i = 0 to n - 1 from the top to the bottom.
 * lower and upper must be padded at lower[0] and upper[n - 1].
 *
 * The vectors are modified in-place. The solution is returned to rhs.
 */
inline void SolveTridiagonalSystem(std::vector<double>* lower,
                                   std::vector<double>* diag,
                                   std::vector<double>* upper,
                                   std::vector<double>* rhs) {
  CHECK_NOTNULL(lower);
  CHECK_NOTNULL(diag);
  CHECK_NOTNULL(upper);
  CHECK_NOTNULL(rhs);

  const std::size_t dim = rhs->size();

  CHECK_NE(dim, 0);
  CHECK_EQ(lower->size(), dim);
  CHECK_EQ(diag->size(), dim);
  CHECK_EQ(upper->size(), dim);

  // Elimination step.
  for (std::size_t i = 0; i < dim - 1; ++i) {
    if (std::abs(diag->at(i)) >= std::abs(lower->at(i + 1))) {
      // Normal tridiagonal algorithm.
      const double w = lower->at(i + 1) / diag->at(i);
      diag->at(i + 1) -= w * upper->at(i);
      rhs->at(i + 1) -= w * rhs->at(i);
      lower->at(i + 1) = 0.0;
    } else {
      // Pivoting. Here, we interchange the i-th row and the (i+1)-th
      // row, then eliminate. Unlike the other branch, the lower
      // triangular element lower[i+1] will remain. This affects the
      // back substitution step below.
      const double w = diag->at(i) / lower->at(i + 1);

      const double u = diag->at(i + 1);
      diag->at(i) = lower->at(i + 1);
      diag->at(i + 1) = upper->at(i) - w * u;
      lower->at(i + 1) = upper->at(i + 1);
      upper->at(i + 1) *= -w;
      upper->at(i) = u;

      const double r = rhs->at(i);
      rhs->at(i) = rhs->at(i + 1);
      rhs->at(i + 1) = r - w * rhs->at(i + 1);
    }
  }

  // Back-substitution step.
  rhs->at(dim - 1) /= diag->at(dim - 1);

  for (std::size_t i_rev = 2; i_rev <= dim; ++i_rev) {
    const std::size_t i = dim - i_rev;
    if (i == dim - 2) {
      rhs->at(i) -= upper->at(i) * rhs->at(i + 1);
    } else {
      rhs->at(i) -=
          upper->at(i) * rhs->at(i + 1) + lower->at(i + 1) * rhs->at(i + 2);
    }
    rhs->at(i) /= diag->at(i);
  }
}

/*
 * Cubic spline functor for interpolating one-dimensional series of values.
 * The original code of can be found in https://github.com/snsinfu/cxx-spline/.
 */
class CubicSpline {
  // Order of the polynomial.
  static constexpr int kOrder = 3;

  // A SplineData defines a single piece of spline function:
  //   f(t) = sum( a_k * (t - t_0)^k , 0 <= k <= 3 ).
  // t_0 is the knot and a_0,...,a_3 are the coefficients.
  struct SplineData {
    double knot;
    std::array<double, kOrder + 1> coefficients;
  };

 public:
  /*
   * Specifies the boundary conditions used to determine the splines.
   */
  enum class BoundaryCondition {
    kNatural = 0,
    kNotAKnot = 1,
  };

  /*
   * Constructs a cubic spline function that passes through given knots.
   */
  CubicSpline(const std::vector<double>& t, const std::vector<double>& x,
              BoundaryCondition bc = BoundaryCondition::kNatural) {
    MakeSpline(t, x, bc);
    MakeBins();
  }

  /*
   * Evaluates the cubic splines at `t`.
   */
  double operator()(double t) const {
    const SplineData& spline = FindSpline(t);

    double value = spline.coefficients[kOrder];
    for (int i = kOrder - 1; i >= 0; --i) {
      value *= t - spline.knot;
      value += spline.coefficients[i];
    }

    return value;
  }

  /*
   * Evaluates the derivatives of cubic splines at `t`.
   */
  double deriv(int order, double t) const {
    const SplineData& spline = FindSpline(t);

    double value = 0.0;
    switch (order) {
      case 1: {
        value = 3.0 * spline.coefficients[kOrder];
        for (int i = kOrder - 1; i >= 1; --i) {
          value *= t - spline.knot;
          value += i * spline.coefficients[i];
        }
        break;
      }
      case 2: {
        value = 6.0 * spline.coefficients[kOrder] * (t - spline.knot) +
                2.0 * spline.coefficients[kOrder - 1];
        break;
      }
      case 3: {
        value = 6.0 * spline.coefficients[kOrder];
        break;
      }
      default: {
        break;
      }
    }
    return value;
  }

 private:
  /*
   * Constructs spline segments `splines_` for given set of knot points.
   */
  void MakeSpline(const std::vector<double>& knots,
                  const std::vector<double>& values, BoundaryCondition bc) {
    if (knots.size() != values.size()) {
      throw std::invalid_argument("input lengths mismatch");
    }
    if (knots.size() <= 1) {
      throw std::invalid_argument("insufficient number of knots");
    }

    const std::size_t n_knots = knots.size();
    const std::size_t n_segments = n_knots - 1;

    std::vector<double> intervals(n_segments);
    std::vector<double> slopes(n_segments);

    for (std::size_t i = 0; i < n_segments; ++i) {
      const double dt = knots[i + 1] - knots[i];
      const double dx = values[i + 1] - values[i];

      if (dt <= 0.0) {
        throw std::invalid_argument("knots must be strictly increasing");
      }

      intervals[i] = dt;
      slopes[i] = dx / dt;
    }

    // Let M[i] be the second derivative of the i-th spline at the i-th
    // knot, i = 0,...,n where n is the number of segments. The vector M is
    // given by a tridiagonal system of equations:
    //
    // [ D[0] U[0]                      ] [ M[0]   ]   [ Y[0]   ]
    // [ L[1] D[1] U[1]                 ] [ M[1]   ]   [ Y[1]   ]
    // [      L[2] D[2] U[2]            ] [ M[2]   ]   [ Y[2]   ]
    // [      ...  ...  ...             ] [   :    ] = [   :    ]
    // [           L[n-1] D[n-1] U[n-1] ] [ M[n-1] ]   [ Y[n-1] ]
    // [                  L[n]   D[n]   ] [ M[n]   ]   [ Y[n]   ]
    //
    // We define the coefficients L, D, U and Y below and solve the
    // equations for M.

    std::vector<double> L(n_segments + 1);
    std::vector<double> D(n_segments + 1);
    std::vector<double> U(n_segments + 1);
    std::vector<double> Y(n_segments + 1);

    if (n_segments == 1) {
      // Natural (which gives a straight line) is the only sensible choice
      // if there is only one spline.
      bc = BoundaryCondition::kNatural;
    }

    switch (bc) {
      case BoundaryCondition::kNatural:
        // These coefficients are derived by assuming that the boundary
        // splines are linear.
        D[0] = 1.0;
        D[n_segments] = 1.0;
        break;

      case BoundaryCondition::kNotAKnot:
        // These coefficients are derived by assuming that the boundary
        // splines are identical to their adjacent splines.
        {
          const double h0 = intervals[0];
          const double h1 = intervals[1];
          const double s0 = slopes[0];
          const double s1 = slopes[1];
          D[0] = h0 - h1;
          U[0] = 2.0 * h0 + h1;
          Y[0] = -6.0 * h0 / (h0 + h1) * (s0 - s1);
        }
        {
          const double h0 = intervals[n_segments - 1];
          const double h1 = intervals[n_segments - 2];
          const double s0 = slopes[n_segments - 1];
          const double s1 = slopes[n_segments - 2];
          D[n_segments] = h0 - h1;
          L[n_segments] = 2.0 * h0 + h1;
          Y[n_segments] = 6.0 * h0 / (h0 + h1) * (s0 - s1);
        }
        break;
    }

    // The remaining coefficients are derived from the spline conditions.
    for (std::size_t i = 1; i < n_segments; ++i) {
      L[i] = intervals[i - 1];
      D[i] = 2.0 * (intervals[i - 1] + intervals[i]);
      U[i] = intervals[i];
      Y[i] = 6.0 * (slopes[i] - slopes[i - 1]);
    }

    // Solve the equations. Solution is returned to Y.
    SolveTridiagonalSystem(&L, &D, &U, &Y);
    const std::vector<double>& M = Y;

    // Derive the polynomial coefficients of each spline segment from the
    // second derivatives `M`.
    splines_.clear();
    splines_.reserve(n_segments);

    for (std::size_t i = 0; i < n_segments; ++i) {
      SplineData spline;
      spline.knot = knots[i];
      spline.coefficients[0] = values[i];
      spline.coefficients[1] =
          slopes[i] - (M[i + 1] + 2.0 * M[i]) * intervals[i] / 6.0;
      spline.coefficients[2] = M[i] / 2.0;
      spline.coefficients[3] = (M[i + 1] - M[i]) / (6.0 * intervals[i]);
      splines_.push_back(spline);
    }

    splines_.shrink_to_fit();

    lower_bound_ = knots.front();
    upper_bound_ = knots.back();
  }

  /*
   * Builds a bin-based index `bins_` that is used to quickly find a spline
   * segment covering a query point.
   */
  void MakeBins() {
    bin_interval_ =
        (upper_bound_ - lower_bound_) / static_cast<double>(splines_.size());

    // We need to map uniformly-spaced bins to spline segments that may
    // not be uniformly spaced. Here `index` identifies a spline segment.
    // We iterate over bins and synchronize spline segments accordingly.
    std::size_t index = 0;

    for (int bin = 0;; ++bin) {
      const double bin_edge = lower_bound_ + bin_interval_ * bin;

      while (index + 1 < splines_.size() &&
             splines_[index + 1].knot <= bin_edge) {
        ++index;
      }

      // Now, the lower edge of the bin is covered by the spline segment
      // at `index`.
      bins_.push_back(index);

      if (index + 1 == splines_.size()) {
        break;
      }
    }

    bins_.shrink_to_fit();
  }

  /*
   * Returns the spline segment that covers given point `t`.
   */
  const SplineData& FindSpline(double t) const {
    if (t <= lower_bound_) {
      return splines_.front();
    }
    if (t >= upper_bound_) {
      return splines_.back();
    }

    auto bin = static_cast<std::size_t>((t - lower_bound_) / bin_interval_);
    if (bin >= bins_.size()) {
      bin = bins_.size() - 1;
    }
    std::size_t index = bins_[bin];

    CHECK_GE(t, splines_[index].knot);

    for (; index + 1 < splines_.size(); ++index) {
      if (t < splines_[index + 1].knot) {
        break;
      }
    }

    return splines_[index];
  }

  std::vector<SplineData> splines_;
  std::vector<std::size_t> bins_;
  double lower_bound_;
  double upper_bound_;
  double bin_interval_;
};

}  // namespace stanley_controller
