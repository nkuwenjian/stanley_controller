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

// band matrix solver
class band_matrix {
 private:
  std::vector<std::vector<double>> m_upper;  // upper band
  std::vector<std::vector<double>> m_lower;  // lower band
 public:
  band_matrix() = delete;                  // constructor
  band_matrix(int dim, int n_u, int n_l);  // constructor
  virtual ~band_matrix() = default;        // destructor
  void resize(int dim, int n_u, int n_l);  // init with dim,n_u,n_l
  int dim() const;                         // matrix dimension
  int num_upper() const { return static_cast<int>(m_upper.size()) - 1; }
  int num_lower() const { return static_cast<int>(m_lower.size()) - 1; }
  // access operator
  double& operator()(int i, int j);       // write
  double operator()(int i, int j) const;  // read
  // we can store an additional diagonal (in m_lower)
  double& saved_diag(int i);
  double saved_diag(int i) const;
  void lu_decompose();
  std::vector<double> r_solve(const std::vector<double>& b) const;
  std::vector<double> l_solve(const std::vector<double>& b) const;
  std::vector<double> lu_solve(const std::vector<double>& b,
                               bool is_lu_decomposed = false);
};

}  // namespace stanley_controller
