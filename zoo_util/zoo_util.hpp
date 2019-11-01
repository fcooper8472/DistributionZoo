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

#ifndef ZOO_UTIL_HPP_
#define ZOO_UTIL_HPP_

#include <numeric>
#include <tuple>

namespace zoo {

template <class real> std::tuple<real, real> moments(const std::vector<real> &sample) {

  const real mean = std::accumulate(sample.begin(), sample.end(), real{0.0}) / sample.size();
  const real var =
      std::inner_product(sample.begin(), sample.end(), sample.begin(), real{0.0}) / sample.size() -
      mean * mean;

  return std::make_tuple(mean, var);
}

template <class real> real median(std::vector<real> &sample) {

  const auto half_way = sample.size() / 2;
  std::nth_element(sample.begin(), std::next(sample.begin(), half_way), sample.end());
  return sample.at(half_way);
}

} // namespace zoo

#endif // ZOO_UTIL_HPP_
