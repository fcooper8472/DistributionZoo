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

#include "zoo_util.hpp"

#define REAL_TYPES float, double, long double

TEMPLATE_TEST_CASE("Utility tests", "[util]", REAL_TYPES) {

  const TestType e = std::numeric_limits<TestType>::epsilon() * 1000;

  // Vector is just full of 1s
  std::vector<TestType> v1(11, TestType{1.0});
  CHECK(v1.size() == 11ul);

  const auto [mean1, var1] = zoo::moments(v1);
  CHECK(mean1 == Approx(TestType{1.0}).epsilon(e));
  CHECK(var1 == Approx(TestType{0.0}).epsilon(e));

  const auto median1 = zoo::median(v1);
  CHECK(median1 == Approx(TestType{1.0}).epsilon(e));

  // Vector has nontrivial contents
  const auto x1 = TestType{1.23};
  const auto x2 = TestType{2.34};
  const auto x3 = TestType{7.89};
  const auto x4 = TestType{-2.67};
  const auto x5 = TestType{-13.4};

  std::vector<TestType> v2 = {x1, x2, x3, x4, x5};
  CHECK(v2.size() == 5ul);

  const auto hand_mean = (x1 + x2 + x3 + x4 + x5) / TestType{5.0};
  const auto hand_var =
      (x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5) / TestType{5.0} - hand_mean * hand_mean;
  const auto hand_median = TestType{1.23};

  const auto [mean2, var2] = zoo::moments(v2);
  CHECK(mean2 == Approx(hand_mean).epsilon(e));
  CHECK(var2 == Approx(hand_var).epsilon(e));

  const auto median2 = zoo::median(v2);
  CHECK(median2 == Approx(hand_median).epsilon(e));
}
