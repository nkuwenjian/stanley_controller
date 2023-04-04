/******************************************************************************
 * Copyright (c) 2022, NKU Mobile & Flying Robotics Lab
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

#include "stanley_controller/common/cubic_spline_interpolation.h"

#include <chrono>  // NOLINT

#include "glog/logging.h"

#include "stanley_controller/third_party/spline/spline.h"

namespace stanley_controller {

void CubicSplineInterpolation::Interpolate(const std::vector<double>& knot_x,
                                           const std::vector<double>& knot_y,
                                           std::vector<double>* spline_x,
                                           std::vector<double>* spline_y,
                                           std::vector<double>* spline_yaw) {
  const auto start_timestamp = std::chrono::system_clock::now();

  // setup auxiliary "time grid"
  double tmin = 0.0;
  double tmax = 0.0;
  std::vector<double> T;

  CreateTimeGrid(&T, &tmin, &tmax, knot_x, knot_y);

  // define a spline for each coordinate x, y
  tk::spline sx;
  tk::spline sy;
  sx.set_points(T, knot_x);
  sy.set_points(T, knot_y);

  // evaluates spline and outputs data to be used with gnuplot
  const double kInterpolationResolution = 0.1;
  const int step =
      static_cast<int>(std::floor((tmax - tmin) / kInterpolationResolution));
  spline_x->resize(step + 1);
  spline_y->resize(step + 1);
  spline_yaw->resize(step + 1);
  for (int i = 0; i <= step; ++i) {
    const double t = tmin + i * kInterpolationResolution;
    spline_x->at(i) = sx(t);
    spline_y->at(i) = sy(t);
    const double deriv_x = sx.deriv(1, t);
    const double deriv_y = sy.deriv(1, t);
    spline_yaw->at(i) = atan2(deriv_y, deriv_x);
  }

  const auto end_timestamp = std::chrono::system_clock::now();
  const std::chrono::duration<double> diff = end_timestamp - start_timestamp;
  LOG(INFO) << "Cubic spline interpolation costs: " << diff.count() * 1e3
            << " ms";
}

void CubicSplineInterpolation::CreateTimeGrid(std::vector<double>* T,
                                              double* tmin, double* tmax,
                                              const std::vector<double>& X,
                                              const std::vector<double>& Y) {
  CHECK_EQ(X.size(), Y.size());
  CHECK_GT(X.size(), 2);

  // setup a "time variable" so that we can interpolate x and y
  // coordinates as a function of time: (X(t), Y(t))
  T->resize(X.size());
  T->front() = 0.0;
  for (size_t i = 1; i < T->size(); i++) {
    // time is proportional to the distance, i.e. we go at a const speed
    T->at(i) = T->at(i - 1) + hypot(X[i] - X[i - 1], Y[i] - Y[i - 1]);
  }

  *tmin = T->front();
  *tmax = T->back();
}

}  // namespace stanley_controller
