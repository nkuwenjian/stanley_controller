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

#include <fstream>
#include <string>
#include <vector>

#include "glog/logging.h"

#include "stanley_controller/common/cubic_spline_interpolation.h"
#include "stanley_controller/common/state.h"

using stanley_controller::CubicSplineInterpolation;
using stanley_controller::kControlPeriod;
using stanley_controller::StanleyController;
using stanley_controller::State;

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = true;
  CHECK_EQ(argc, 1);

  std::vector<double> x{0.0, 100.0, 100.0, 50.0, 60.0};
  std::vector<double> y{0.0, 0.0, -30.0, -20.0, 0.0};

  std::vector<double> spline_x;
  std::vector<double> spline_y;
  std::vector<double> spline_yaw;
  CubicSplineInterpolation::Interpolate(x, y, &spline_x, &spline_y,
                                        &spline_yaw);
  CHECK_EQ(spline_x.size(), spline_y.size());
  CHECK_EQ(spline_x.size(), spline_yaw.size());
  LOG(INFO) << "size of cubic spline: " << spline_x.size();

  // Initial state
  State state(0.0, 5.0, M_PI / 9.0, 0.0);
  const size_t last_idx = spline_x.size() - 1;
  size_t target_idx = 0;
  StanleyController::CalcTargetIndex(state, spline_x, spline_y, &target_idx,
                                     nullptr);
  const double max_simulation_time = 100.0;
  const double target_speed = 30.0 / 3.6;  // [m/s]
  double time = 0.0;
  std::vector<double> t{0.0};
  std::vector<double> result_x{state.x()};
  std::vector<double> result_y{state.y()};
  std::vector<double> result_v{state.v()};

  while (max_simulation_time >= time && last_idx > target_idx) {
    double acc = StanleyController::PidControl(state.v(), target_speed);
    double steer = 0.0;
    size_t current_target_idx = 0;
    StanleyController::StanleyControl(state, spline_x, spline_y, spline_yaw,
                                      target_idx, &steer, &current_target_idx);
    state.Update(acc, steer);
    time += kControlPeriod;
    target_idx = current_target_idx;
    t.push_back(time);
    result_x.push_back(state.x());
    result_y.push_back(state.y());
    result_v.push_back(state.v());
  }
  CHECK_GE(last_idx, target_idx) << "Cannot reach goal";

  // Log spline data for debug
  std::string file = "../data/spline.txt";
  std::ofstream fout;
  fout.open(file);
  if (!fout.is_open()) {
    LOG(ERROR) << "Failed to open " << file;
    return 1;
  }
  for (size_t i = 0; i < spline_x.size(); ++i) {
    fout << std::fixed << spline_x[i] << "," << spline_y[i] << ","
         << spline_yaw[i] << "\n";
  }
  fout.close();

  // Log tracking result for debug
  file = "../data/track.txt";
  fout.open(file);
  if (!fout.is_open()) {
    LOG(ERROR) << "Failed to open " << file;
    return 1;
  }
  for (size_t i = 0; i < result_x.size(); ++i) {
    fout << std::fixed << t[i] << "," << result_x[i] << "," << result_y[i]
         << "," << result_v[i] << "\n";
  }
  fout.close();

  return 0;
}
