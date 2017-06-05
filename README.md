# Mean shift clustering

A simple C++ implementation of the mean shift clustering algorithm with the following properties and restrictions:

- Comprised of a single header file for the main algorithm, with additional headers for extra utilities.
- Can use OpenMP to parallelize the execution.
- Templated to support different scalar types.
- Supports custom kernel, distance and bandwidth estimator functions.
- Template specialization system to support any container/structure type.
- Requires C++11 support.

## Basic usage

The most important function is `mean_shift_cluster`. The following is a simple example illustrating its usage:

```cpp
List<Point<Scalar>> points = get_points();
std::vector<msc::Cluster<Scalar>> clusters = msc::mean_shift_cluster<Scalar>(
    std::begin(points), std::end(points), 3, metric, kernel, estimator);
```

The main algorithm is in `msh.h`. Implementations of some common metrics and kernels and a test estimator is in `msc.metrics.h`, `msc.kernels.h` and `msc.estimators.h`, respectively, but any custom functor or function that implements the appropriate signature is valid. The signatures are:

```cpp
double metric(const Scalar* a, const Scalar* b, int dim);
double kernel(double distance);
double estimator(const Scalar* point,
    InputIterator first, InputIterator last, int dim, Metric metric);
```

To account for custom containers for point data, an `msc::Accessor` functor can be specialized to the desired type and the compiler will pick the correct specialization to access the container values. The signature is:

```cpp
struct Accessor<Scalar, Point>
{
    static const Scalar* data(const Point& point);
}
```

Specializations for common types are included in `msc.accessors.h`, so that things like raw arrays, `std::array`s and `std::vector`s work just by including this header (Note: it must be included after the main `msc.h` header for the compiler to know about the base template beforehand).

A helper header `msc` can be used to include all these headers in a single line.

## Tests and examples

A generic calculator is included in the file `main.cpp` that reads points from a file or the standard input and dumps the clustered points to the standard output. Some tests are included in the following files:

- `test_custom_struct`: Exemplifies the use of a custom structure (`Point3`) to store points.
- `test_1d_flat_vector`: Uses a flat vector to store 1D points. This configuration works thanks to one of the accessors included in `msc.accessors.h`.

The script `test.sh` uses the main executable to read the dataset in `test.txt` (obtained from [here](http://www.uni-marburg.de/fb12/arbeitsgruppen/datenbionik/data)) and plots the results with gnuplot.
