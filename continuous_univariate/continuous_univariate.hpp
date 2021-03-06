/*
MIT License

Copyright (c) 2019 University of Oxford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef CONTINUOUS_UNIVARIATE_HPP_
#define CONTINUOUS_UNIVARIATE_HPP_

#include <cmath>
#include <random>
#include <vector>

namespace zoo {

template <class real> constexpr real pi = real{3.14159265358979323846264338L};

template <class real> class ContinuousUnivariate {
private:
  std::random_device mRd{};

protected:
  std::mt19937 mMt{mRd()};

public:
  virtual real pdf(real x) = 0;
  virtual real log_pdf(real x) = 0;
  virtual real rand() = 0;
  virtual std::vector<real> randn(std::size_t n) {
    std::vector<real> sample(n);
    for (auto &x : sample) {
      x = this->rand();
    }
    return sample;
  }
};

template <class real> class Beta : public ContinuousUnivariate<real> {
private:
  // Params
  real mAlpha;
  real mBeta;

  // Dists
  std::gamma_distribution<real> mDistX;
  std::gamma_distribution<real> mDistY;

  // Cached constants for Pdf & LogPdf
  real m1OnBetaFn;
  real mLogBetaFn;
  real mAm1;
  real mBm1;

public:
  explicit Beta(const real alpha = 1.0, const real beta = 1.0) : mAlpha(alpha), mBeta(beta) {

    // Both params must be positive
    assert(mAlpha > real{0.0});
    assert(mBeta > real{0.0});

    mDistX = std::gamma_distribution<real>{mAlpha, real{1.0}};
    mDistY = std::gamma_distribution<real>{mBeta, real{1.0}};

    // Constants for Beta function evaluations
    m1OnBetaFn = std::tgamma(mAlpha + mBeta) / (std::tgamma(mAlpha) * std::tgamma(mBeta));
    mLogBetaFn = std::lgamma(mAlpha + mBeta) - (std::lgamma(mAlpha) + std::lgamma(mBeta));

    // Other useful constants
    mAm1 = mAlpha - real{1.0};
    mBm1 = mBeta - real{1.0};
  }

  real pdf(const real x) override {
    if (x > real{0.0} && x < real{1.0}) {
      return std::pow(x, mAm1) * std::pow(real{1.0} - x, mBm1) * m1OnBetaFn;
    } else {
      return real{0.0};
    }
  }

  real log_pdf(const real x) override {
    if (x > real{0.0} && x < real{1.0}) {
      return mAm1 * std::log(x) + mBm1 * std::log1p(-x) + mLogBetaFn;
    } else {
      return -std::numeric_limits<real>::infinity();
    }
  }

  real rand() override {
    const real x = mDistX(this->mMt);
    const real y = mDistY(this->mMt);
    return x / (x + y);
  }
};

template <class real> class Normal : public ContinuousUnivariate<real> {
private:
  // Params
  real mMean;
  real mStdDev;

  // Dist
  std::normal_distribution<real> mDist;

  // Cached constants for Pdf & LogPdf
  real m2SigSq;
  real mPrefactor;
  real mLogPrefactor;

public:
  explicit Normal(const real mean = 0.0, const real std_dev = 1.0) : mMean(mean), mStdDev(std_dev) {

    // Standard deviation must be positive
    assert(mStdDev > real{0.0});

    mDist = std::normal_distribution<real>{mMean, mStdDev};

    m2SigSq = real{2.0} * mStdDev * mStdDev;
    mPrefactor = real{1.0} / std::sqrt(zoo::pi<real> * m2SigSq);
    mLogPrefactor = real{-0.5} * std::log(zoo::pi<real> * m2SigSq);
  }

  real pdf(const real x) override {
    return mPrefactor * std::exp(-(x - mMean) * (x - mMean) / m2SigSq);
  }

  real log_pdf(const real x) override {
    return mLogPrefactor - (x - mMean) * (x - mMean) / m2SigSq;
  }

  real rand() override { return mDist(this->mMt); }
};

} // namespace zoo

#endif // CONTINUOUS_UNIVARIATE_HPP_
