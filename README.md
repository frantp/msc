# Mean shift clustering

A simple C++ implementation of the mean shift clustering algorithm with the following properties and restrictions:

- Comprised of a single header file.
- Can use OpenMP to parallelize the execution.
- Templated to support different scalar types.
- Supports custom kernel and distance functions.
- Requires C++11 support.
- Expects points as `std::vector`s of the scalar type.

A simple test is included in the file `main.cpp` that reads points in CSV format from the standard input and dumps the clustered points to the standard output. The script `test.sh` uses the main executable to read the dataset in `test.csv` and plots the results with gnuplot.
