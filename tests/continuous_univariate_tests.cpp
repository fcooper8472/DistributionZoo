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

#include "catch.hpp"

#include <iostream>
#include <limits>

#include "continuous_univariate.hpp"
#include "zoo_util.hpp"

#define REAL_TYPES float, double, long double

TEMPLATE_TEST_CASE("Beta values", "[beta]", REAL_TYPES) {

  const TestType e = std::numeric_limits<TestType>::epsilon() * 1000;

  const TestType alpha = 2.6L;
  const TestType beta = 4.9L;
  zoo::Beta<TestType> dist{alpha, beta};

  // Regular PDF
  CHECK(dist.pdf(-1.0) == Approx(TestType{0.0L}).epsilon(e));
  CHECK(dist.pdf(0.5) == Approx(TestType{1.399459344806713569240L}).epsilon(e));
  CHECK(dist.pdf(2.0) == Approx(TestType{0.0L}).epsilon(e));

  // Log PDF
  CHECK(std::isinf(dist.log_pdf(-1.0)));
  CHECK(dist.log_pdf(0.5) == Approx(TestType{0.3360859797527134507530L}).epsilon(e));
  CHECK(std::isinf(dist.log_pdf(2.0)));

  // Sample
  const TestType big_e{0.1};
  const std::size_t n = 10001;
  auto sample = dist.randn(n);

  const auto [mean, var] = zoo::moments(sample);
  CHECK(mean == Approx(alpha / (alpha + beta)).epsilon(big_e));

  const auto hand_var =
      alpha * beta / ((alpha + beta) * (alpha + beta) * (alpha + beta + TestType{1.0}));
  CHECK(var == Approx(hand_var).epsilon(big_e));

  const auto median = zoo::median(sample);
  const auto hand_median = (alpha - TestType{1.0 / 3.0}) / (alpha + beta - TestType{2.0 / 3.0});
  CHECK(median == Approx(hand_median).epsilon(big_e));
}

TEMPLATE_TEST_CASE("Normal values", "[normal]", REAL_TYPES) {

  const TestType e = std::numeric_limits<TestType>::epsilon() * 1000;

  const TestType dist_mean{8.9L};
  const TestType dist_std_dev{2.3L};
  zoo::Normal<TestType> dist{dist_mean, dist_std_dev};

  // Regular PDF
  CHECK(dist.pdf(5.0) == Approx(TestType{0.04119387068037555522332L}).epsilon(e));
  CHECK(dist.pdf(9.6) == Approx(TestType{0.1656030770867795793562L}).epsilon(e));

  // Log PDF
  CHECK(dist.log_pdf(5.0) == Approx(TestType{-3.189465803587791871442L}).epsilon(e));
  CHECK(dist.log_pdf(9.6) == Approx(TestType{-1.798161455761704914921L}).epsilon(e));

  // Sample
  const TestType big_e{0.1};
  const std::size_t n = 10001;
  auto sample = dist.randn(n);

  const auto [mean, var] = zoo::moments(sample);
  CHECK(mean == Approx(dist_mean).epsilon(big_e));
  CHECK(var == Approx(dist_std_dev * dist_std_dev).epsilon(big_e));

  const auto median = zoo::median(sample);
  CHECK(median == Approx(dist_mean).epsilon(big_e));
}
