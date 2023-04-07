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

#pragma once

#include <vector>

namespace stanley_controller {

class CubicSpline {
 public:
  // spline types
  enum spline_type {
    linear = 10,          // linear interpolation
    cspline = 30,         // cubic splines (classical C^2)
    cspline_hermite = 31  // cubic hermite splines (local, only C^1)
  };

  // boundary condition type for the spline end-points
  enum bd_type { first_deriv = 1, second_deriv = 2, not_a_knot = 3 };

  CubicSpline() = default;
  CubicSpline(const std::vector<double>& X, const std::vector<double>& Y,
              spline_type type = cspline, bd_type left = second_deriv,
              double left_value = 0.0, bd_type right = second_deriv,
              double right_value = 0.0);
  virtual ~CubicSpline() = default;

  void set_points(const std::vector<double>& x, const std::vector<double>& y,
                  spline_type type = cspline);

  // evaluates the spline at point x
  double operator()(double x) const;
  double deriv(int order, double x) const;

 private:
  std::vector<double> m_x;  // x coordinates of points
  std::vector<double> m_y;  // y coordinates of points
  // interpolation parameters
  // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
  // where a_i = y_i, or else it won't go through grid points
  std::vector<double> m_b, m_c, m_d;  // spline coefficients
  double m_c0;                        // for left extrapolation
  spline_type m_type = cspline;
  bd_type m_left = second_deriv;
  bd_type m_right = second_deriv;
  double m_left_value = 0.0;
  double m_right_value = 0.0;

  void set_coeffs_from_b();             // calculate c_i, d_i from b_i
  size_t find_closest(double x) const;  // closest idx so that m_x[idx]<=x
};

}  // namespace stanley_controller
