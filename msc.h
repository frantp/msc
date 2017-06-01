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
#include <iterator>
#include <type_traits>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace msc
{
template <class T>
struct Cluster
{
    std::vector<T> mode;
    std::vector<std::size_t> members;

    inline Cluster(const T* mode, int dim)
        : mode(mode, mode + dim), members() {}
};

template <class T, class C>
struct Accessor
{
    inline static const T* data(const C& point)
    {
        static_assert(False<T>::value, "Accessor not implemented for container");
    }

private:
    template <class>
    struct False : std::false_type {};
};

template <class T, class ForwardIterator,
          class Metric, class Kernel, class Estimator>
inline std::vector<T> meanshift(const T* point,
    ForwardIterator first, ForwardIterator last, int dim,
    Metric metric, Kernel kernel, Estimator estimator)
{
    if (dim <= 0)
        throw std::invalid_argument("Dimension must be greater than 0");
    typedef typename std::iterator_traits<ForwardIterator>::value_type C;
    std::vector<T> shifted(dim);
    const auto ibw = estimator(point, first, last, dim, metric);
    double total_weight = 0;

    for (auto it = first; it != last; it++)
    {
        const T* pt = Accessor<T, C>::data(*it);
        const auto dist = metric(pt, point, dim);
        const auto weight = kernel(dist * ibw);
        for (int k = 0; k < shifted.size(); k++)
            shifted[k] += pt[k] * weight;
        total_weight += weight;
    }

    for (int k = 0; k < shifted.size(); k++)
        shifted[k] /= total_weight;

    return shifted;
}

template <class T, class ForwardIterator,
          class Metric, class Kernel, class Estimator>
inline std::vector<std::vector<T>> meanshift(
    ForwardIterator first, ForwardIterator last, int dim,
    Metric metric, Kernel kernel, Estimator estimator,
    double epsilon = std::numeric_limits<float>::epsilon(),
    int max_iter = std::numeric_limits<int>::max())
{
    if (dim <= 0)
        throw std::invalid_argument("Dimension must be greater than 0");
    typedef typename std::iterator_traits<ForwardIterator>::value_type C;
    std::vector<std::vector<T>> shifted(std::distance(first, last));
    std::size_t i = 0;
    for (auto it = first; it != last; it++, i++)
    {
        const T* pt = Accessor<T, C>::data(*it);
        shifted[i].resize(dim);
        for (int k = 0; k < dim; k++)
            shifted[i][k] = pt[k];
    }

    #pragma omp parallel for
    for (std::size_t i = 0; i < shifted.size(); i++)
    {
        const T* pt = shifted[i].data();
        int iter = 0;
        double d = 0;
        do
        {
            const auto point = meanshift(
                pt, first, last, dim, metric, kernel, estimator);
            d = metric(pt, point.data(), dim);
            shifted[i] = point;
            iter++;
        }
        while (d > epsilon && iter < max_iter);
    }

    return shifted;
}

template <class T, class InputIterator, class Metric>
inline std::vector<Cluster<T>> cluster(
    InputIterator first, InputIterator last, int dim, Metric metric,
    double epsilon = std::numeric_limits<float>::epsilon())
{
    if (dim <= 0)
        throw std::invalid_argument("Dimension must be greater than 0");
    typedef typename std::iterator_traits<InputIterator>::value_type C;
    std::vector<Cluster<T>> clusters;
    std::size_t i = 0;
    for (auto it = first; it != last; it++, i++)
    {
        const T* pt = Accessor<T, C>::data(*it);
        int c = 0;
        for (; c < clusters.size(); c++)
            if (metric(pt, clusters[c].mode.data(), dim) <= epsilon)
                break;
        if (c == clusters.size())
            clusters.emplace_back(pt, dim);
        clusters[c].members.emplace_back(i);
    }

    return clusters;
}

template <class T, class ForwardIterator,
          class Metric, class Kernel, class Estimator>
inline std::vector<Cluster<T>> meanshiftcluster(
    ForwardIterator first, ForwardIterator last, int dim,
    Metric metric, Kernel kernel, Estimator estimator,
    double epsilon = std::numeric_limits<float>::epsilon(),
    int max_iter = std::numeric_limits<int>::max())
{
    const auto shifted = meanshift<T>(
        first, last, dim, metric, kernel, estimator, epsilon, max_iter);
    return cluster<T>(
        std::begin(shifted), std::end(shifted), dim, metric, epsilon);
}
} // namespace msc
