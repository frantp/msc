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

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <iostream>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "meanshift.h"

typedef double Scalar;

std::vector<std::vector<Scalar>> load(std::istream& in);
void dump(const std::vector<std::vector<Scalar>>& points,
    const std::vector<msc::Cluster<Scalar>>& clusters);

int main(int argc, char** argv)
{
    std::istream* in = &std::cin;
    std::ifstream infile;
    #ifdef _WIN32
    if (_isatty(_fileno(stdin)))
    #else
    if (isatty(fileno(stdin)))
    #endif
    {
        std::cout << "Input CSV file: ";
        std::string filename;
        std::cin >> filename;
        infile.open(filename);
        in = &infile;
    }

    const double bandwidth = argc > 1 ? std::stof(argv[1]) : 0.1;
    std::cerr << "Kernel bandwidth: " << bandwidth << std::endl;
    const auto points = load(*in);
    std::cerr << "Num. points: " << points.size() << std::endl;
    const auto clusters = msc::cluster(points, msc::kernel::Gaussian(bandwidth));
    std::cerr << "Clusters:" << std::endl;
    for (const auto& cluster : clusters)
    {
        std::cerr << " - Num. elems.: " << cluster.members.size() << "; Center:";
        for (const auto& value : cluster.center)
            std::cerr << ", " << value;
        std::cerr << std::endl;
    }
    dump(points, clusters);
    return 0;
}

std::vector<std::vector<Scalar>> load(std::istream& in)
{
    std::vector<std::vector<Scalar>> points;
    std::string line;
    while (std::getline(in, line))
    {
        points.emplace_back();
        auto& point = points.back();
        std::istringstream lin(line);
        std::string token;
        while (std::getline(lin, token, ','))
            point.push_back(std::stod(token));
    }
    return points;
}

void dump(const std::vector<std::vector<Scalar>>& points,
    const std::vector<msc::Cluster<Scalar>>& clusters)
{
    for (std::size_t c = 0; c < clusters.size(); c++)
    {
        for (const auto& index : clusters[c].members)
        {
            std::cout << c;
            const auto point = points[index];
            for (std::size_t i = 0; i < point.size(); i++)
                std::cout << "," << point[i];
            std::cout << std::endl;
        }
    }
}
