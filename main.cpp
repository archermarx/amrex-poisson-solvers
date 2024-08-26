#include <AMReX.H>
#include <AMReX_Print.H>

int main (int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
#if 0
  {
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
  }
#endif 
  amrex::Finalize();
}
