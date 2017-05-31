# Mean shift clustering

A simple C++ implementation of the mean shift clustering algorithm with the following properties and restrictions:

- Comprised of a single header file.
- Can use OpenMP to parallelize the execution.
- Templated to support different scalar types.
- Supports custom kernel, distance and bandwidth estimator functions.
- Template specialization system to support any container/structure type.
- Requires C++11 support.

The main algorithm is in `msh.h`. Implementations of some common metrics and kernels and a test estimator is in `msc.metrics.h`, `msc.kernels.h` and `msc.estimators.h`, respectively, but any custom functor or function that implements the appropriate signature is valid. The signatures are:

```cpp
double metric(const T* a, const T* b, int dim);
double kernel(double d);
double estimator(const T* p, InputIterator first, InputIterator last, int dim, Metric metric);
```

To account for custom containers for point data, an `msc::Accessor` functor can be specialized to the desired type and the algorithm will pick the correct specialization when it needs to access the container values. The signature is:

```cpp
struct Accessor<Scalar, Container>
{
    static Scalar* operator()(const Container& container);
}
```

Specializations for common types are included in `msc.accessors.h`, so that raw arrays, `std::array`s and `std::vector`s work just by including this header (Note: it must be included after the main `msc.h` header for the compiler to know about the base template beforehand).

A simple test is included in the file `main.cpp` that reads points in CSV format from a file or the standard input and dumps the clustered points to the standard output. It exemplifies the use of a custom structure (`Point3`) to store points. The script `test.sh` uses the main executable to read the dataset in `test.csv` and plots the results with gnuplot.
