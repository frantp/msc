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

#include <array>
#include <vector>

namespace msc
{
template <class T>
struct Accessor<T, T>
{
    inline static const T* data(const T& point)
    {
        return &point;
    }
};

template <class T>
struct Accessor<T, T*>
{
    inline static const T* data(const T* point)
    {
        return point;
    }
};

template <class T, int N>
struct Accessor<T, std::array<T, N>>
{
    inline static const T* data(const std::array<T, N>& point)
    {
        return point.data();
    }
};

template <class T, class Alloc>
struct Accessor<T, std::vector<T, Alloc>>
{
    inline static const T* data(const std::vector<T, Alloc>& point)
    {
        return point.data();
    }
};
} // namespace msc
