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

#include "msc"

#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <iostream>
#include <chrono>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

template <class T>
struct Point3
{
    T x, y, z;
};

namespace msc
{
template <class T>
struct Accessor<T, Point3<T>>
{
    inline static const T* data(const Point3<T>& container)
    {
        return &container.x;
    }
};
} // namespace msc

typedef double Scalar;
typedef Point3<Scalar> Container;

std::vector<Container> load(std::istream& in);
void dump(const std::vector<Container>& points,
    const std::vector<msc::Cluster<Scalar>>& clusters);

int main()
{
    double bandwidth = 3;
    std::ifstream in("test.txt");
    std::cerr << "Kernel bandwidth: " << bandwidth << std::endl;
    const auto points = load(in);
    std::cerr << "Num. points: " << points.size() << std::endl;
    if (points.empty())
        return 0;
    const auto t0 = std::chrono::high_resolution_clock::now();
    const auto clusters = msc::mean_shift_cluster<Scalar>(
        std::begin(points), std::end(points),
        sizeof(Container) / sizeof(Scalar),
        msc::metrics::L2Sq(),
        msc::kernels::ParabolicSq(),
        msc::estimators::Constant(bandwidth));
    const auto t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "Clusters (" << clusters.size() << "):" << std::endl;
    for (const auto& cluster : clusters)
    {
        std::cerr << " - Num. elems.: " << cluster.members.size() << "; Mode:";
        for (const auto& value : cluster.mode)
            std::cerr << " " << value;
        std::cerr << std::endl;
    }
    std::cerr << "Elapsed time: " << std::chrono::duration_cast<
        std::chrono::microseconds>(t1 - t0).count() / 1e6 << std::endl;
    dump(points, clusters);
    return 0;
}

std::vector<Container> load(std::istream& in)
{
    std::vector<Container> points;
    std::string line;
    while (std::getline(in, line))
    {
        points.emplace_back();
        auto& point = points.back();
        std::istringstream lin(line);
        std::string token;
        lin >> token;
        point.x = std::stod(token);
        lin >> token;
        point.y = std::stod(token);
        lin >> token;
        point.z = std::stod(token);
    }
    return points;
}

void dump(const std::vector<Container>& points,
    const std::vector<msc::Cluster<Scalar>>& clusters)
{
    for (std::size_t c = 0; c < clusters.size(); c++)
    {
        for (const auto& index : clusters[c].members)
        {
            const auto& point = points[index];
            std::cout << c << " "
                      << point.x << " "
                      << point.y << " "
                      << point.z << std::endl;
        }
    }
}
