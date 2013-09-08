A library for the approximation of the elastodynamics equation 
==============================================================

Installation on Linux
---------------------
The compilation of the library requires CMake. To install, run

    git clone https://github.com/carlomr/tspeed.git  
    cd tspeed  
    mkdir build  
    cmake ..  
    make
    make install  
    make doc  

This will by default install the library in `/usr/local/lib` and the header files in `/usr/local/include/tspeed`. To choose a different installation directory, use

    cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/dir

Two test are generated as `Examples/Lamb` and `Examples/Wedge`. For the source generating them, see `Examples/src/Wedge.cpp` and `Examples/src/Lamb.cpp`.




