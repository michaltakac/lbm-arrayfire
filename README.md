Lattice-Boltzmann code leveraging ArrayFire
=====

## Building The Project

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

### ArrayFire Support and Contact Info

* Google Groups: https://groups.google.com/forum/#!forum/arrayfire-users
* ArrayFire Services:  [Consulting](http://arrayfire.com/consulting/)  |  [Support](http://arrayfire.com/support/)   |  [Training](http://arrayfire.com/training/)
* ArrayFire Blogs: http://arrayfire.com/blog/
* Email: <mailto:technical@arrayfire.com>
