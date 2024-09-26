cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=RelWithDbgInfo \
  -DAMReX_OMP=ON \
  -DAMReX_MPI=OFF \
  -DAMReX_SPACEDIM=2

cmake --build build -j 8

