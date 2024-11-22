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

#include "stanley_controller/common/spline.h"

#include <iostream>
#include <vector>

int main() {
  {
    // Some points (t,y) on a curve y = f(t)
    std::vector<double> t = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 2.0, 1.0, 2.0};

    // Spline interpolation (and extrapolation) of the points
    stanley_controller::CubicSpline spline(t, y);

    std::cout << spline(0.5) << "\n";  // 1.44976
    std::cout << spline(1.5) << "\n";  // 2.65072
    std::cout << spline(6.0) << "\n";  // 3
  }
  {
    std::vector<double> X = {0.1, 0.4, 1.2, 1.8, 2.0};  // must be increasing
    std::vector<double> Y = {0.1, 0.7, 0.6, 1.1, 0.9};

    stanley_controller::CubicSpline spline(X, Y);
    double t = 1.5;
    std::cout << spline(t) << "\n";                         // 0.915345
    std::cout << std::fixed << spline.deriv(1, t) << "\n";  // 1.223601
  }
  {
    std::vector<double> X = {0.1, 0.4};
    std::vector<double> Y = {0.1, 0.7};

    stanley_controller::CubicSpline spline(X, Y);
    std::cout << spline(0.25) << "\n";  // 0.4
  }

  return 0;
}
