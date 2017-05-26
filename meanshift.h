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
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace msc
{
template <typename Scalar>
struct Cluster
{
    std::vector<Scalar> center;
    std::vector<std::size_t> members;

    inline Cluster(const std::vector<Scalar>& center)
        : center(center), members() {}
};

namespace kernel
{
struct Gaussian
{
    inline Gaussian(double bandwidth)
        : ibw2_(1.0 / (bandwidth * bandwidth)) {}

    inline double operator()(double d) const
    {
        return std::exp(-0.5 * d * d * ibw2_);
    }

private:
    double ibw2_;
};
} // namespace kernel

namespace dist
{
struct SqEuclidean
{
    template <typename Scalar>
    inline double operator()(
        const std::vector<Scalar>& a,
        const std::vector<Scalar>& b) const
    {
        Scalar d = 0;
        for (int i = 0; i < a.size(); i++)
            d += (a[i] - b[i]) * (a[i] - b[i]);
        return d;
    }
};
} // namespace dist

namespace detail
{
double EPS = std::numeric_limits<float>::epsilon();

template <typename Scalar, typename Kernel, typename Distance>
inline std::vector<Scalar> meanshift(const std::vector<Scalar>& point,
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Distance distance)
{
    std::vector<Scalar> shifted(point.size());
    double total_weight = 0;

    #pragma omp parallel for
    for (std::size_t p = 0; p < points.size(); p++)
    {
        const auto weight = kernel(distance(point, points[p]));
        for (int i = 0; i < shifted.size(); i++)
            shifted[i] += points[p][i] * weight;
        total_weight += weight;
    }

    for (int i = 0; i < shifted.size(); i++)
        shifted[i] /= total_weight;
    return shifted;
}

template <typename Scalar, typename Distance>
inline std::vector<Cluster<Scalar>> cluster(
    const std::vector<std::vector<Scalar>>& shifted, Distance distance)
{
    std::vector<Cluster<Scalar>> clusters;

    for (int i = 0; i < shifted.size(); i++)
    {
        int c = 0;
        for (; c < clusters.size(); c++)
            if (distance(shifted[i], clusters[c].center) <= EPS)
                break;
        if (c == clusters.size())
            clusters.emplace_back(shifted[i]);
        clusters[c].members.emplace_back(i);
    }

    return clusters;
}
} // namespace detail

template <typename Scalar, typename Kernel, typename Distance>
inline std::vector<std::vector<Scalar>> meanshift(
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Distance distance)
{
    std::vector<std::vector<Scalar>> shifted(points);

    #pragma omp parallel for
    for (std::size_t i = 0; i < shifted.size(); i++)
    {
        double d = 0;
        do
        {
            const auto p = detail::meanshift(
                shifted[i], points, kernel, distance);
            d = distance(p, shifted[i]);
            shifted[i] = p;
        }
        while (d > detail::EPS);
    }

    return shifted;
}

template <typename Scalar, typename Kernel,
    typename Distance = dist::SqEuclidean>
inline std::vector<Cluster<Scalar>> cluster(
    const std::vector<std::vector<Scalar>>& points,
    Kernel kernel, Distance distance = Distance())
{
    return detail::cluster(meanshift(points, kernel, distance), distance);
}

} // namespace msc
