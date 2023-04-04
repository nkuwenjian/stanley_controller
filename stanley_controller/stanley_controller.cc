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

#include "stanley_controller/stanley_controller.h"

#include <algorithm>
#include <limits>

namespace stanley_controller {

void StanleyController::CalcTargetIndex(const State& state,
                                        const std::vector<double>& ref_x,
                                        const std::vector<double>& ref_y,
                                        size_t* target_idx,
                                        double* error_front_axle) {
  CHECK_NOTNULL(target_idx);
  CHECK_EQ(ref_x.size(), ref_y.size());
  CHECK(!ref_x.empty());

  // Calc front axle position
  const double front_x = state.x() + kWheelBase * cos(state.yaw());
  const double front_y = state.y() + kWheelBase * sin(state.yaw());

  // Search nearest point index
  *target_idx = 0;
  double min_dist = std::numeric_limits<double>::max();
  for (size_t i = 0; i < ref_x.size(); ++i) {
    const double dx = front_x - ref_x[i];
    const double dy = front_y - ref_y[i];
    const double dxy = hypot(dx, dy);
    if (dxy < min_dist) {
      min_dist = dxy;
      *target_idx = i;
    }
  }

  // Project RMS error onto front axle vector
  const double front_axle_vec_x = -cos(state.yaw() + M_PI / 2.0);
  const double front_axle_vec_y = -sin(state.yaw() + M_PI / 2.0);
  if (error_front_axle != nullptr) {
    *error_front_axle = (front_x - ref_x[*target_idx]) * front_axle_vec_x +
                        (front_y - ref_y[*target_idx]) * front_axle_vec_y;
  }
}

void StanleyController::StanleyControl(const State& state,
                                       const std::vector<double>& ref_x,
                                       const std::vector<double>& ref_y,
                                       const std::vector<double>& ref_yaw,
                                       size_t last_target_idx, double* steer,
                                       size_t* current_target_idx) {
  double error_front_axle;
  CalcTargetIndex(state, ref_x, ref_y, current_target_idx, &error_front_axle);
  *current_target_idx = std::max(last_target_idx, *current_target_idx);

  // theta_e corrects the heading error
  const double theta_e =
      NormalizeAngle(ref_yaw[*current_target_idx] - state.yaw());
  // theta_d corrects the cross track error
  constexpr double kStanleyControllerGain = 0.5;
  const double theta_d =
      atan2(kStanleyControllerGain * error_front_axle, state.v());
  // Steering control
  *steer = Clamp(theta_e + theta_d, -kMaxSteer, kMaxSteer);
}

}  // namespace stanley_controller
