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

#include <vector>
#include <cmath>

namespace msc
{
namespace metrics
{
struct L1
{
    template <typename Scalar>
    inline double operator()(const Scalar* a, const Scalar* b, int dim) const
    {
        Scalar d = 0;
        for (int i = 0; i < dim; i++)
            d += a[i] - b[i];
        return d;
    }
};

struct L2
{
    template <typename Scalar>
    inline double operator()(const Scalar* a, const Scalar* b, int dim) const
    {
        Scalar d = 0;
        for (int i = 0; i < dim; i++)
            d += (a[i] - b[i]) * (a[i] - b[i]);
        return std::sqrt(d);
    }
};

struct L2Sq
{
    template <typename Scalar>
    inline double operator()(const Scalar* a, const Scalar* b, int dim) const
    {
        Scalar d = 0;
        for (int i = 0; i < dim; i++)
            d += (a[i] - b[i]) * (a[i] - b[i]);
        return d;
    }
};

struct Inf
{
    template <typename Scalar>
    inline double operator()(const Scalar* a, const Scalar* b, int dim) const
    {
        Scalar d = 0;
        for (int i = 0; i < dim; i++)
        {
            const auto t = a[i] - b[i];
            if (d < t) d = t;
        }
        return d;
    }
};
} // namespace metrics
} // namespace msc
