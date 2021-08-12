Lattice-Boltzmann code leveraging ArrayFire
=====

Paper: https://doi.org/10.3390/math9151793

If you use the code from this repository, please cite this paper:

```
@article{takac2021lbm,
    author = {Takáč, Michal and Petráš, Ivo},
    title = {Cross-Platform GPU-Based Implementation of Lattice Boltzmann Method Solver Using ArrayFire Library},
    journal = {Mathematics},
    volume = {9},
    year = {2021},
    number = {15},
    article-number = {1793},
    url = {https://www.mdpi.com/2227-7390/9/15/1793},
    issn = {2227-7390},
    abstract = {This paper deals with the design and implementation of cross-platform, D2Q9-BGK and D3Q27-MRT, lattice Boltzmann method solver for 2D and 3D flows developed with ArrayFire library for high-performance computing. The solver leverages ArrayFire’s just-in-time compilation engine for compiling high-level code into optimized kernels for both CUDA and OpenCL GPU backends. We also provide C++ and Rust implementations and show that it is possible to produce fast cross-platform lattice Boltzmann method simulations with minimal code, effectively less than 90 lines of code. An illustrative benchmarks (lid-driven cavity and Kármán vortex street) for single and double precision floating-point simulations on 4 different GPUs are provided.},
    doi = {10.3390/math9151793}
}
```

## Building C++ code

Once ArrayFire is installed on your machine, compile the project
using `cmake` and `make`:

    mkdir build
    cd build
    cmake ..
    make

If ArrayFire is not installed to a system directory, you will need to specify
the directory which contains the `ArrayFireConfig.cmake` as an argument to the
`cmake` invocation. This configuration file is located within the
`share/ArrayFire` subdirectory of the ArrayFire installation. For example,
if you were to install ArrayFire to the `local` directory within your home
folder, the invocation of `cmake` above would be replaced with the following:

    cmake -DArrayFire_DIR=$HOME/local/share/ArrayFire/cmake ..

Note to self: `cmake -DArrayFire_DIR=/opt/arrayfire/share/ArrayFire/cmake ..`

## Building Rust code

Run the examples by running

    cd rust
    cargo run --example channel

or

    cd rust
    cargo run --example lid

### ArrayFire Support and Contact Info

* Google Groups: https://groups.google.com/forum/#!forum/arrayfire-users
* ArrayFire Services:  [Consulting](http://arrayfire.com/consulting/)  |  [Support](http://arrayfire.com/support/)   |  [Training](http://arrayfire.com/training/)
* ArrayFire Blogs: http://arrayfire.com/blog/
* Email: <mailto:technical@arrayfire.com>
