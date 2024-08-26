cmake -S . -B build \
  -DAMReX_OMP=ON \
  -DAMReX_MPI=OFF \
  -DCMAKE_BUILD_TYPE=RelWithDbgInfo

cmake --build build -j 8

