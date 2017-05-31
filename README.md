# Mean shift clustering

A simple C++ implementation of the mean shift clustering algorithm with the following properties and restrictions:

- Comprised of a single header file.
- Can use OpenMP to parallelize the execution.
- Templated to support different scalar types.
- Supports custom kernel, distance and bandwidth estimator functions.
- The Point container/structure must implement `operator[]` and the data must be contiguous in memory. This means that, for example, c arrays, `std::array`s and `std::vector`s are all valid containers.
- Requires C++11 support.

The main algorithm is in `msh.h`. Implementations of some common metrics, kernels and a test estimator is in `msc.metrics.h`, `msc.kernels.h` and `msc.estimators.h`, respectively, but any custom functor or function that implements the appropriate signature is valid. The signatures are:

```
double metric(const T* a, const T* b, int dim);
double kernel(double d);
double estimator(const T* p, InputIterator first, InputIterator last, int dim, Metric metric);
```

A simple test is included in the file `main.cpp` that reads points in CSV format from a file or the standard input and dumps the clustered points to the standard output. The script `test.sh` uses the main executable to read the dataset in `test.csv` and plots the results with gnuplot.
