#include <AMReX.H>
#include "PoissonProblem.H"

int main (int argc, char *argv[]) {
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");
        PoissonProblem myprob;
        myprob.solve();
        myprob.compute_norms();
        myprob.write_plotfile();
    }
    amrex::Finalize();
}
