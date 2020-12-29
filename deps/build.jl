using libcxxwrap_julia_jll

# assume you have correctly install libint2

cxxwrapdir = libcxxwrap_julia_jll.artifact_dir
run(`clang++ -O3 -std=c++17 -I "$cxxwrapdir/include" -I "/usr/include/julia" coulomb.cpp -shared -fPIC -o libcoulomb.so -ljulia -lint2 -pthread -L"$cxxwrapdir/lib" -lcxxwrap_julia`)