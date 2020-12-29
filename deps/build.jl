using CxxWrap

cxxwrapdir = CxxWrap.prefix_path()
run(`g++ -O3 -std=c++17 -I "$cxxwrapdir/include" -I "/usr/include/julia" -I"/usr/include/eigen3" coulomb.cpp -shared -fPIC -o libcoulomb.so -L"$cxxwrapdir/lib" -lcxxwrap_julia`)