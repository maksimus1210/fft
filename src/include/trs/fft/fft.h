#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <complex>
#include <assert.h>
namespace fft {

typedef std::complex<double> complex_double;
constexpr double PI = M_PI;
constexpr double DOUBLE_PI = 2 * PI;

template <typename It> void fft(It f, It l) {
  typedef size_t size_type;
  auto n = l - f;
  assert((n - 1) & n != 0);
  std::vector<size_type> rev(n);
  size_type log_n = std::log2(n);
  size_type log_n_minus_one = log_n - 1;

  #pragma omp parallel for
  for (size_type i = 0; i < n; ++i) {
    for (size_type j = 0; j < log_n; ++j) {
      if (i & (size_type(1) << j)) {
        rev[i] |= size_type(1) << size_type(log_n_minus_one - j);
      }
    }
  }

  std::vector<complex_double> wlen_pw(n);

  size_type i = 0;
  auto it = f;
  for (; it != l;++it, ++i) {
    if (i < rev[i]) {
      iter_swap(f + rev[i], it);
    }
  }

    wlen_pw[0] = complex_double(1.0, 0.0);
  for (size_type len = 2; len <= n; len <<= 1) {
    double ang = DOUBLE_PI / len ;
    auto half = len >> 1;

    complex_double wlen(cos(ang), sin(ang));
    for (size_type i = 1; i < half; ++i) {
      wlen_pw[i] = wlen_pw[i - 1] * wlen;
    }

    #pragma omp parallel for 
    for (size_type i = 0; i < n; i += len) {
      complex_double t;
      auto pu = f + i;
      auto pv = f + i + half;
      auto pu_end = f + i + half;
      auto pw = std::begin(wlen_pw);
      for (; pu != pu_end; ++pu, ++pv, ++pw) {
        t = *pv * *pw;
        *pv = *pu - t;
        *pu += t;
      }
    }
  }

}

template <typename It> void fft_inv(It f, It l) {
  // typedef typename std::iterator_traits<It>::size_type size_type;
  typedef size_t size_type;
  auto n = std::distance(f, l);
  assert((n - 1) & n != 0);
  std::vector<size_type> rev(n);
  size_type log_n = std::log2(n);
  size_type log_n_minus_one = log_n - 1;
  #pragma omp parallel for 
  for (size_type i = 0; i < n; ++i) {
    for (size_type j = 0; j < log_n; ++j) {
      if (i & (size_type(1) << j)) {
        rev[i] |= size_type(1) << size_type(log_n_minus_one - j);
      }
    }
  }

  std::vector<complex_double> wlen_pw(n, complex_double(1.0, 0.0));
  size_type i = 0;
  auto it = f;
  while (it != l) {
    if (i < rev[i]) {
      iter_swap(f + rev[i], it);
    }
    ++it;
    ++i;
  }
  it = f;
  //complex_double one_zero(1.0, 0.0);
  for (size_type len = 2; len <= n; len <<= 1) {
    double ang = -2.0 * PI / len;
    auto half = len >> 1;
    complex_double wlen(cos(ang), sin(ang));
    for (size_type i = 1; i < half; ++i) {
      wlen_pw[i] = wlen_pw[i - 1] * wlen;
    }

    #pragma omp parallel for 
    for (size_type i = 0; i < n; i += len) {
      complex_double t;
      auto pu = f + i;
      auto pv = f + i + half;
      auto pu_end = f + i + half;
      auto pw = std::begin(wlen_pw);
      for (; pu != pu_end; ++pu, ++pv, ++pw) {
        t = *pv * *pw;
        *pv = *pu - t;
        *pu += t;
      }
    }
  }

  #pragma omp parallel for 
  for (size_type i = 0; i < n; i++) { *(f + i) /= n; };
}

}
