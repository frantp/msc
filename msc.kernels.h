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
namespace kernels
{
struct Uniform
{
    inline double operator()(double d) const
    {
        return d <= 1 ? 1 : 0;
    }
};

struct Triangular
{
    inline double operator()(double d) const
    {
        return d <= 1 ? 1 - std::abs(d) : 0;
    }
};

struct Parabolic
{
    inline double operator()(double d) const
    {
        return d <= 1 ? 1 - d * d : 0;
    }
};

struct ParabolicSq
{
    inline double operator()(double d2) const
    {
        return d2 <= 1 ? 1 - d2 : 0;
    }
};

struct Biweight
{
    inline double operator()(double d) const
    {
        const auto x = 1 - d * d;
        return d <= 1 ? x * x : 0;
    }
};

struct BiweightSq
{
    inline double operator()(double d2) const
    {
        const auto x = 1 - d2;
        return d2 <= 1 ? x * x : 0;
    }
};

struct Triweight
{
    inline double operator()(double d) const
    {
        const auto x = 1 - d * d;
        return d <= 1 ? x * x * x : 0;
    }
};

struct TriweightSq
{
    inline double operator()(double d2) const
    {
        const auto x = 1 - d2;
        return d2 <= 1 ? x * x * x : 0;
    }
};

struct Tricube
{
    inline double operator()(double d) const
    {
        const auto x = 1 - d * d * d;
        return d <= 1 ? x * x * x : 0;
    }
};

struct TricubeCu
{
    inline double operator()(double d3) const
    {
        const auto x = 1 - d3;
        return d3 <= 1 ? x * x * x : 0;
    }
};

struct Gaussian
{
    inline double operator()(double d) const
    {
        return std::exp(-0.5 * d * d);
    }
};

struct GaussianSq
{
    inline double operator()(double d2) const
    {
        return std::exp(-0.5 * d2);
    }
};

struct Cosine
{
    inline double operator()(double d) const
    {
        return d <= 1 ? std::cos(M_PI_2 * d) : 0;
    }
};

struct Logistic
{
    inline double operator()(double d) const
    {
        return 1.0 / (2 + std::exp(d) + std::exp(-d));
    }
};

struct Sigmoid
{
    inline double operator()(double d) const
    {
        return 1.0 / (std::exp(d) + std::exp(-d));
    }
};

struct Silverman
{
    inline double operator()(double d) const
    {
        const auto x = M_SQRT1_2 * std::abs(d);
        return std::exp(-x) * std::sin(x + M_PI_4);
    }
};
} // namespace kernels
} // namespace msc
