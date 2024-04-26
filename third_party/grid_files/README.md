# C++ Grid File Library

A C++ single-file header-only implementation of in-memory multi-dimensional Grid Files. It doesn't require any  third party libraries. Just include the **GridFile.h** file in your code.

Grid Files are a data structure that is particularly useful as a multidimensional index. They offer quick data access and range queries through a multidimensional indexing scheme. Their memory overhead is much less than unlike other multidimensional structures such as R-Trees or Kd-Trees.

### Background
Due to the scarcity of open-source implementations of multidimensional indexes that support range queries in C++, this library is to serve as a baseline for range query benchmarks. 


### Library Features
- bulk loading capability of very large datasets in a column-layout fashion. (see examples)
- range query operation with a visitor pattern.
- run-time dynamic definition of dataset dimensionality
- per grid cell column-layout definition instead of row-layout to reduce pointers memory-overhead.

### Lack of Multidimensional Indexes as Baselines
- The range query functionality supported by this library is also offered in Boost RTree library.
- It's useful to test the performance of the two libraries. While R-Tree may suffer from memory overhead for very large datasets, this implementation of grid files can serve as an alternative.

- It is unfortunate that many recent scientific works in the field of multidimensional indexes use closed source baselines when testing the performances of their advanced indexes. 
- Hence, this work serves as an attempt to increase the open-source baseline implementations.
