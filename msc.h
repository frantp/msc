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
#include <memory>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace msc
{
constexpr double EPS = std::numeric_limits<float>::epsilon();

template <class T, class Alloc = std::allocator<T>>
struct Cluster
{
    std::vector<T, Alloc> mode;
    std::vector<std::size_t, Alloc> members;

    inline Cluster(const std::vector<T>& mode)
        : mode(mode), members() {}
};

template <class T, class InputIterator,
          class Kernel, class Metric, class Estimator>
inline std::vector<T> meanshift(const std::vector<T>& point,
    InputIterator first, InputIterator last,
    Kernel kernel, Metric metric, Estimator estimator)
{
    std::vector<T> shifted(point.size());
    const auto ibw = estimator(point, first, last, metric);
    double total_weight = 0;

    for (auto it = first; it != last; it++)
    {
        const auto& pt = *it;
        const auto dist = metric(point.data(), pt.data(), point.size());
        const auto weight = kernel(dist * ibw);
        for (int i = 0; i < shifted.size(); i++)
            shifted[i] += pt[i] * weight;
        total_weight += weight;
    }

    for (int i = 0; i < shifted.size(); i++)
        shifted[i] /= total_weight;
    return shifted;
}

template <class T, class ForwardIterator,
          class Kernel, class Metric, class Estimator>
inline std::vector<std::vector<T>> meanshift(
    ForwardIterator first, ForwardIterator last,
    Kernel kernel, Metric metric, Estimator estimator)
{
    std::vector<std::vector<T>> shifted(first, last);

    #pragma omp parallel for
    for (std::size_t p = 0; p < shifted.size(); p++)
    {
        double d = 0;
        do
        {
            const auto point = meanshift(
                shifted[p], first, last, kernel, metric, estimator);
            d = metric(point.data(), shifted[p].data(), point.size());
            shifted[p] = point;
        }
        while (d > EPS);
    }

    return shifted;
}

template <class T, class Alloc = std::allocator<T>, class InputIterator,
          class Metric>
inline std::vector<Cluster<T, Alloc>> cluster(
    InputIterator first, InputIterator last, Metric metric)
{
    std::vector<Cluster<T, Alloc>> clusters;
    std::size_t i = 0;
    for (auto it = first; it != last; it++, i++)
    {
        const auto& pt = *it;
        int c = 0;
        for (; c < clusters.size(); c++)
            if (metric(pt.data(), clusters[c].mode.data(), pt.size()) <= EPS)
                break;
        if (c == clusters.size())
            clusters.emplace_back(pt);
        clusters[c].members.emplace_back(i);
    }

    return clusters;
}

template <class T, class Alloc = std::allocator<T>, class ForwardIterator,
          class Kernel, class Metric, class Estimator>
inline std::vector<Cluster<T, Alloc>> meanshiftcluster(
    ForwardIterator first, ForwardIterator last,
    Kernel kernel, Metric metric, Estimator estimator)
{
    const auto shifted = meanshift<T>(first, last, kernel, metric, estimator);
    return cluster<T>(std::begin(shifted), std::end(shifted), metric);
}
} // namespace msc
