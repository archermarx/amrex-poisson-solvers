cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=RelWithDbgInfo \
  -DAMReX_MPI=OFF \
  -DAMReX_SPACEDIM=3 \
  -DAMReX_TINY_PROFILE=ON \
  -DAMReX_BASE_PROFILE=OFF

cmake --build build -j 20

