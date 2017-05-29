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
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace msc
{
template <typename Scalar>
struct Cluster
{
    std::vector<Scalar> mode;
    std::vector<std::size_t> members;

    inline Cluster(const std::vector<Scalar>& mode)
        : mode(mode), members() {}
};

namespace detail
{
double EPS = std::numeric_limits<float>::epsilon();

template <typename Scalar, typename Kernel, typename Metric>
inline std::vector<Scalar> meanshift(const std::vector<Scalar>& point,
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Metric metric, Scalar ibw)
{
    std::vector<Scalar> shifted(point.size());
    double total_weight = 0;

    #pragma omp parallel for reduction(+: total_weight)
    for (std::size_t p = 0; p < points.size(); p++)
    {
        const auto dist = metric(point.data(), points[p].data(), point.size());
        const auto weight = kernel(dist * ibw);
        for (int i = 0; i < shifted.size(); i++)
            shifted[i] += points[p][i] * weight;
        total_weight += weight;
    }

    for (int i = 0; i < shifted.size(); i++)
        shifted[i] /= total_weight;
    return shifted;
}

template <typename Scalar, typename Metric>
inline std::vector<Cluster<Scalar>> cluster(
    const std::vector<std::vector<Scalar>>& shifted, Metric metric)
{
    std::vector<Cluster<Scalar>> clusters;

    for (int i = 0; i < shifted.size(); i++)
    {
        int c = 0;
        for (; c < clusters.size(); c++)
            if (metric(shifted[i].data(), clusters[c].mode.data(),
                shifted[i].size()) <= EPS)
                break;
        if (c == clusters.size())
            clusters.emplace_back(shifted[i]);
        clusters[c].members.emplace_back(i);
    }

    return clusters;
}
} // namespace detail

template <typename Scalar, typename Kernel, typename Metric>
inline std::vector<std::vector<Scalar>> meanshift(
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Metric metric, Scalar bandwidth)
{
    std::vector<std::vector<Scalar>> shifted(points);
    const auto ibw = 1 / bandwidth;

    #pragma omp parallel for
    for (std::size_t i = 0; i < shifted.size(); i++)
    {
        double d = 0;
        do
        {
            const auto p = detail::meanshift(
                shifted[i], points, kernel, metric, ibw);
            d = metric(p.data(), shifted[i].data(), p.size());
            shifted[i] = p;
        }
        while (d > detail::EPS);
    }

    return shifted;
}

template <typename Scalar, typename Kernel, typename Metric>
inline std::vector<Cluster<Scalar>> cluster(
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Metric metric, Scalar bandwidth)
{
    return detail::cluster(
        meanshift(points, kernel, metric, bandwidth), metric);
}
} // namespace msc
