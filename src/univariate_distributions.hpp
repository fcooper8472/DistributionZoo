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
class UnivariateDistribution {
public:
  virtual RealType pdf(RealType x) = 0;
  virtual RealType log_pdf(RealType x) = 0;
};

template <typename RealType>
class Normal : public UnivariateDistribution<RealType> {
private:
  // Params
  RealType mMean;
  RealType mStdDev;

  // Cached constants for Pdf & LogPdf
  RealType m2SigSq;
  RealType mPrefactor;
  RealType mLogPrefactor;

public:
  explicit Normal(RealType mean = 0.0, RealType std_dev = 1.0)
      : mMean(mean), mStdDev(std_dev) {

    // Standard deviation must be positive
    assert(mStdDev > 0.0);

    m2SigSq = 2.0 * mStdDev * mStdDev;
    mPrefactor = 1.0 / std::sqrt(M_PI * m2SigSq);
    mLogPrefactor = -0.5 * std::log(M_PI * m2SigSq);
  }

  RealType pdf(const RealType x) override {
    return mPrefactor * std::exp(-(x - mMean) * (x - mMean) / m2SigSq);
  }

  RealType log_pdf(const RealType x) override {
    return mLogPrefactor - (x - mMean) * (x - mMean) / m2SigSq;
  }
};


} // namespace zoo

#endif // UNIVARIATE_DISTRIBUTIONS_HPP_
