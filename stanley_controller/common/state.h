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

#pragma once

#include <cmath>

#include "glog/logging.h"

#include "stanley_controller/common/util.h"

namespace stanley_controller {

static constexpr double kMaxSteer = M_PI / 6.0;  // [rad] max steering angle
static constexpr double kControlPeriod = 0.1;    // [s] Control period
static const double kWheelBase = 2.9;            // [m] Wheel base of vehicle

class State {
 public:
  State(const double x, const double y, const double yaw, const double v)
      : x_(x), y_(y), yaw_(yaw), v_(v) {}
  virtual ~State() = default;

  void Update(const double acc, const double steer) {
    CHECK_GE(steer, -kMaxSteer);
    CHECK_LE(steer, kMaxSteer);

    x_ += v_ * kControlPeriod * cos(yaw_);
    y_ += v_ * kControlPeriod * sin(yaw_);
    yaw_ += v_ * kControlPeriod / kWheelBase * tan(steer);
    yaw_ = NormalizeAngle(yaw_);
    v_ += acc * kControlPeriod;
  }

  double x() const { return x_; }
  double y() const { return y_; }
  double yaw() const { return yaw_; }
  double v() const { return v_; }

 private:
  double x_ = 0.0;
  double y_ = 0.0;
  double yaw_ = 0.0;
  double v_ = 0.0;
};

}  // namespace stanley_controller
