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

// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <limits>

#include "univariate_distributions.hpp"

#define REAL_TYPES float, double, long double

TEMPLATE_TEST_CASE("Beta values", "[beta]", REAL_TYPES) {

  const TestType e = std::numeric_limits<TestType>::epsilon() * 1000;

  const TestType mean = 8.9L;
  const TestType std_dev = 2.3L;
  zoo::Normal<TestType> dist{mean, std_dev};

  // Regular PDF
  CHECK(dist.pdf(5.0) == Approx(static_cast<TestType>(0.04119387068037555522332L)).epsilon(e));
  CHECK(dist.pdf(9.6) == Approx(static_cast<TestType>(0.1656030770867795793562L)).epsilon(e));

  // Log PDF
  CHECK(dist.log_pdf(5.0) == Approx(static_cast<TestType>(-3.189465803587791871442L)).epsilon(e));
  CHECK(dist.log_pdf(9.6) == Approx(static_cast<TestType>(-1.798161455761704914921L)).epsilon(e));
}

}
