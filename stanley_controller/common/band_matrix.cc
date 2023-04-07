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

#include "stanley_controller/common/band_matrix.h"

#include <algorithm>

#include "glog/logging.h"

namespace stanley_controller {

band_matrix::band_matrix(int dim, int n_u, int n_l) { resize(dim, n_u, n_l); }

void band_matrix::resize(int dim, int n_u, int n_l) {
  CHECK_GT(dim, 0);
  CHECK_GE(n_u, 0);
  CHECK_GE(n_l, 0);
  m_upper.resize(n_u + 1);
  m_lower.resize(n_l + 1);
  for (std::vector<double>& upper : m_upper) {
    upper.resize(dim);
  }
  for (std::vector<double>& lower : m_lower) {
    lower.resize(dim);
  }
}

int band_matrix::dim() const {
  return m_upper.empty() ? 0 : static_cast<int>(m_upper[0].size());
}

// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double& band_matrix::operator()(int i, int j) {
  int k = j - i;  // what band is the entry
  CHECK_GE(i, 0);
  CHECK_LT(i, dim());
  CHECK_GE(j, 0);
  CHECK_LT(j, dim());
  CHECK_LE(-num_lower(), k);
  CHECK_LE(k, num_upper());
  // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
  if (k >= 0) {
    return m_upper[k][i];
  }
  return m_lower[-k][i];
}

double band_matrix::operator()(int i, int j) const {
  int k = j - i;  // what band is the entry
  CHECK_GE(i, 0);
  CHECK_LT(i, dim());
  CHECK_GE(j, 0);
  CHECK_LT(j, dim());
  CHECK_LE(-num_lower(), k);
  CHECK_LE(k, num_upper());
  // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
  if (k >= 0) {
    return m_upper[k][i];
  }
  return m_lower[-k][i];
}

// second diag (used in LU decomposition), saved in m_lower
double band_matrix::saved_diag(int i) const {
  CHECK_GE(i, 0);
  CHECK_LT(i, dim());
  return m_lower[0][i];
}

double& band_matrix::saved_diag(int i) {
  CHECK_GE(i, 0);
  CHECK_LT(i, dim());
  return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose() {
  int i_max;
  int j_max;
  int j_min;
  double x;

  // preconditioning
  // normalize column i so that a_ii=1
  for (int i = 0; i < this->dim(); i++) {
    CHECK_NE(this->operator()(i, i), 0.0);
    this->saved_diag(i) = 1.0 / this->operator()(i, i);
    j_min = std::max(0, i - this->num_lower());
    j_max = std::min(this->dim() - 1, i + this->num_upper());
    for (int j = j_min; j <= j_max; j++) {
      this->operator()(i, j) *= this->saved_diag(i);
    }
    this->operator()(i, i) = 1.0;  // prevents rounding errors
  }

  // Gauss LR-Decomposition
  for (int k = 0; k < this->dim(); k++) {
    i_max = std::min(this->dim() - 1,
                     k + this->num_lower());  // num_lower not a mistake!
    for (int i = k + 1; i <= i_max; i++) {
      CHECK_NE(this->operator()(k, k), 0.0);
      x = -this->operator()(i, k) / this->operator()(k, k);
      this->operator()(i, k) = -x;  // assembly part of L
      j_max = std::min(this->dim() - 1, k + this->num_upper());
      for (int j = k + 1; j <= j_max; j++) {
        // assembly part of R
        this->operator()(i, j) =
            this->operator()(i, j) + x * this->operator()(k, j);
      }
    }
  }
}

// solves Ly=b
std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const {
  CHECK_EQ(this->dim(), (int)b.size());
  std::vector<double> x(this->dim());
  int j_start;
  double sum;
  for (int i = 0; i < this->dim(); i++) {
    sum = 0;
    j_start = std::max(0, i - this->num_lower());
    for (int j = j_start; j < i; j++) {
      sum += this->operator()(i, j) * x[j];
    }
    x[i] = (b[i] * this->saved_diag(i)) - sum;
  }
  return x;
}

// solves Rx=y
std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const {
  CHECK_EQ(this->dim(), (int)b.size());
  std::vector<double> x(this->dim());
  int j_stop;
  double sum;
  for (int i = this->dim() - 1; i >= 0; i--) {
    sum = 0;
    j_stop = std::min(this->dim() - 1, i + this->num_upper());
    for (int j = i + 1; j <= j_stop; j++) {
      sum += this->operator()(i, j) * x[j];
    }
    x[i] = (b[i] - sum) / this->operator()(i, i);
  }
  return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
                                          bool is_lu_decomposed) {
  CHECK_EQ(this->dim(), (int)b.size());
  std::vector<double> x;
  std::vector<double> y;
  if (!is_lu_decomposed) {
    this->lu_decompose();
  }
  y = this->l_solve(b);
  x = this->r_solve(y);
  return x;
}

}  // namespace stanley_controller
