/******************************************************************************
 * Copyright (c) 2023, NKU Mobile & Flying Robotics Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Jian Wen (nkuwenjian@gmail.com)
 *****************************************************************************/

#include "stanley_controller/common/cubic_spline.h"

#include <algorithm>

#include "glog/logging.h"

#include "stanley_controller/common/band_matrix.h"

namespace stanley_controller {

CubicSpline::CubicSpline(const std::vector<double>& X,
                         const std::vector<double>& Y, spline_type type,
                         bd_type left, double left_value, bd_type right,
                         double right_value)
    : m_type(type),
      m_left(left),
      m_right(right),
      m_left_value(left_value),
      m_right_value(right_value) {
  set_points(X, Y, m_type);
}

void CubicSpline::set_points(const std::vector<double>& x,
                             const std::vector<double>& y, spline_type type) {
  CHECK_EQ(x.size(), y.size());
  CHECK_GE(x.size(), 3);

  m_type = type;
  m_x = x;
  m_y = y;
  int n = static_cast<int>(x.size());
  // check strict monotonicity of input vector x
  for (int i = 0; i < n - 1; i++) {
    CHECK_LT(m_x[i], m_x[i + 1]);
  }

  if (type == linear) {
    // linear interpolation
    m_d.resize(n);
    m_c.resize(n);
    m_b.resize(n);
    for (int i = 0; i < n - 1; i++) {
      m_d[i] = 0.0;
      m_c[i] = 0.0;
      m_b[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
    }
    // ignore boundary conditions, set slope equal to the last segment
    m_b[n - 1] = m_b[n - 2];
    m_c[n - 1] = 0.0;
    m_d[n - 1] = 0.0;
  } else if (type == cspline) {
    // classical cubic splines which are C^2 (twice cont differentiable)
    // this requires solving an equation system

    // setting up the matrix and right hand side of the equation system
    // for the parameters b[]
    int n_upper = (m_left == CubicSpline::not_a_knot) ? 2 : 1;
    int n_lower = (m_right == CubicSpline::not_a_knot) ? 2 : 1;
    band_matrix A(n, n_upper, n_lower);
    std::vector<double> rhs(n);
    for (int i = 1; i < n - 1; i++) {
      A(i, i - 1) = 1.0 / 3.0 * (x[i] - x[i - 1]);
      A(i, i) = 2.0 / 3.0 * (x[i + 1] - x[i - 1]);
      A(i, i + 1) = 1.0 / 3.0 * (x[i + 1] - x[i]);
      rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
               (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    }
    // boundary conditions
    if (m_left == CubicSpline::second_deriv) {
      // 2*c[0] = f''
      A(0, 0) = 2.0;
      A(0, 1) = 0.0;
      rhs[0] = m_left_value;
    } else {
      CHECK(false);
    }
    if (m_right == CubicSpline::second_deriv) {
      // 2*c[n-1] = f''
      A(n - 1, n - 1) = 2.0;
      A(n - 1, n - 2) = 0.0;
      rhs[n - 1] = m_right_value;
    } else {
      CHECK(false);
    }

    // solve the equation system to obtain the parameters c[]
    m_c = A.lu_solve(rhs);

    // calculate parameters b[] and d[] based on c[]
    m_d.resize(n);
    m_b.resize(n);
    for (int i = 0; i < n - 1; i++) {
      m_d[i] = 1.0 / 3.0 * (m_c[i + 1] - m_c[i]) / (x[i + 1] - x[i]);
      m_b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
               1.0 / 3.0 * (2.0 * m_c[i] + m_c[i + 1]) * (x[i + 1] - x[i]);
    }
    // for the right extrapolation coefficients (zero cubic term)
    // f_{n-1}(x) = y_{n-1} + b*(x-x_{n-1}) + c*(x-x_{n-1})^2
    double h = x[n - 1] - x[n - 2];
    // m_c[n-1] is determined by the boundary condition
    m_d[n - 1] = 0.0;
    m_b[n - 1] = 3.0 * m_d[n - 2] * h * h + 2.0 * m_c[n - 2] * h +
                 m_b[n - 2];  // = f'_{n-2}(x_{n-1})
    if (m_right == first_deriv) {
      m_c[n - 1] = 0.0;  // force linear extrapolation
    }
  } else {
    CHECK(false);
  }

  // for left extrapolation coefficients
  m_c0 = (m_left == first_deriv) ? 0.0 : m_c[0];
}

double CubicSpline::operator()(double x) const {
  // polynomial evaluation using Horner's scheme
  // TODO(all): consider more numerically accurate algorithms, e.g.:
  //   - Clenshaw
  //   - Even-Odd method by A.C.R. Newbery
  //   - Compensated Horner Scheme
  size_t n = m_x.size();
  size_t idx = find_closest(x);

  double h = x - m_x[idx];
  double interpol;
  if (x < m_x[0]) {
    // extrapolation to the left
    interpol = (m_c0 * h + m_b[0]) * h + m_y[0];
  } else if (x > m_x[n - 1]) {
    // extrapolation to the right
    interpol = (m_c[n - 1] * h + m_b[n - 1]) * h + m_y[n - 1];
  } else {
    // interpolation
    interpol = ((m_d[idx] * h + m_c[idx]) * h + m_b[idx]) * h + m_y[idx];
  }
  return interpol;
}

double CubicSpline::deriv(int order, double x) const {
  CHECK_GT(order, 0);
  size_t n = m_x.size();
  size_t idx = find_closest(x);

  double h = x - m_x[idx];
  double interpol;
  if (x < m_x[0]) {
    // extrapolation to the left
    switch (order) {
      case 1:
        interpol = 2.0 * m_c0 * h + m_b[0];
        break;
      case 2:
        interpol = 2.0 * m_c0;
        break;
      default:
        interpol = 0.0;
        break;
    }
  } else if (x > m_x[n - 1]) {
    // extrapolation to the right
    switch (order) {
      case 1:
        interpol = 2.0 * m_c[n - 1] * h + m_b[n - 1];
        break;
      case 2:
        interpol = 2.0 * m_c[n - 1];
        break;
      default:
        interpol = 0.0;
        break;
    }
  } else {
    // interpolation
    switch (order) {
      case 1:
        interpol = (3.0 * m_d[idx] * h + 2.0 * m_c[idx]) * h + m_b[idx];
        break;
      case 2:
        interpol = 6.0 * m_d[idx] * h + 2.0 * m_c[idx];
        break;
      case 3:
        interpol = 6.0 * m_d[idx];
        break;
      default:
        interpol = 0.0;
        break;
    }
  }
  return interpol;
}

// return the closest idx so that m_x[idx] <= x (return 0 if x<m_x[0])
size_t CubicSpline::find_closest(double x) const {
  std::vector<double>::const_iterator it;
  it = std::upper_bound(m_x.begin(), m_x.end(), x);  // *it > x
  size_t idx = std::max(static_cast<int>(std::distance(m_x.begin(), it)) - 1,
                        0);  // m_x[idx] <= x
  return idx;
}

}  // namespace stanley_controller
