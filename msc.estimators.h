// Copyright (c) 2017 Francisco Troncoso Pastoriza
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <cmath>

namespace msc
{
namespace estimators
{
struct Constant
{
    inline explicit Constant(double factor)
        : factor_(1 / factor) {}

    template <class T, class InputIterator, class Metric>
    inline double operator()(const T* point,
        InputIterator first, InputIterator last, int dim, Metric metric) const
    {
        return factor_;
    }

private:
    double factor_;
};

struct MinMaxDistance
{
    inline explicit MinMaxDistance(double factor)
        : factor_(1 / (factor * factor)) {}

    template <class T, class InputIterator, class Metric>
    inline double operator()(const T* point,
        InputIterator first, InputIterator last, int dim, Metric metric) const
    {
        double dmin = 0, dmax = 0;
        for (auto it = first; it != last; it++)
        {
            const auto& pt = &(*it)[0];
            const auto d = metric(pt, point, dim);
            if (d > 0)
            {
                if (d < dmin || dmin == 0) dmin = d;
                if (d > dmax) dmax = d;
            }
        }
        return std::pow(dmax / dmin, factor_) * factor_;
    }

private:
    double factor_;
};
} // namespace estimators
} // namespace msc
