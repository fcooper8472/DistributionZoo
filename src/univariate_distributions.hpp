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

#ifndef UNIVARIATE_DISTRIBUTIONS_HPP_
#define UNIVARIATE_DISTRIBUTIONS_HPP_

#include <cmath>

namespace zoo {

template <typename RealType>
constexpr RealType PI() {
	return static_cast<RealType>(3.14159265358979323846264338L);
}

template <typename RealType> class UnivariateDistribution {
public:
  virtual RealType pdf(RealType x) = 0;
  virtual RealType log_pdf(RealType x) = 0;
};

template <typename RealType> class Normal : public UnivariateDistribution<RealType> {
private:
  // Params
  RealType mMean;
  RealType mStdDev;

  // Cached constants for Pdf & LogPdf
  RealType m2SigSq;
  RealType mPrefactor;
  RealType mLogPrefactor;

public:
  explicit Normal(const RealType mean = 0.0, const RealType std_dev = 1.0)
      : mMean(mean), mStdDev(std_dev) {

    // Standard deviation must be positive
    assert(mStdDev > 0.0);

    m2SigSq = 2.0 * mStdDev * mStdDev;
    mPrefactor = 1.0 / std::sqrt(zoo::PI<RealType>() * m2SigSq);
    mLogPrefactor = -0.5 * std::log(zoo::PI<RealType>() * m2SigSq);
  }

  RealType pdf(const RealType x) override {
    return mPrefactor * std::exp(-(x - mMean) * (x - mMean) / m2SigSq);
  }

  RealType log_pdf(const RealType x) override {
    return mLogPrefactor - (x - mMean) * (x - mMean) / m2SigSq;
  }
};

template <typename RealType> class Beta : public UnivariateDistribution<RealType> {
private:
  // Params
  RealType mAlpha;
  RealType mBeta;

  // Cached constants for Pdf & LogPdf
  RealType m1OnBetaFn;
  RealType mLogBetaFn;
  RealType mAm1;
  RealType mBm1;

public:
  explicit Beta(const RealType alpha = 1.0, const RealType beta = 1.0)
      : mAlpha(alpha), mBeta(beta) {

    // Both params must be positive
    assert(mAlpha > 0.0);
    assert(mBeta > 0.0);

    // Constants for Beta function evaluations
    m1OnBetaFn = std::tgamma(mAlpha + mBeta) / (std::tgamma(mAlpha) * std::tgamma(mBeta));
    mLogBetaFn = std::lgamma(mAlpha + mBeta) - (std::lgamma(mAlpha) + std::lgamma(mBeta));

    // Other useful constants
    mAm1 = mAlpha - 1.0;
    mBm1 = mBeta - 1.0;
  }

  RealType pdf(const RealType x) override {
    if (x > 0.0 && x < 1.0) {
      return std::pow(x, mAm1) * std::pow(1.0 - x, mBm1) * m1OnBetaFn;
    } else {
      return 0.0;
    }
  }

  RealType log_pdf(const RealType x) override {
    if (x > 0.0 && x < 1.0) {
      return mAm1 * std::log(x) + mBm1 * std::log1p(-x) + mLogBetaFn;
    } else {
      return -std::numeric_limits<RealType>::infinity();
    }
  }
};

} // namespace zoo

#endif // UNIVARIATE_DISTRIBUTIONS_HPP_
