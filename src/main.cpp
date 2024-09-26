#include <AMReX.H>
#include "PoissonProblem.H"

int main (int argc, char *argv[]) {

    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");
        auto strt_time = amrex::second();
        PoissonProblem myprob;
        myprob.solve();
        myprob.compute_norms();
        myprob.write_plotfile();
        auto stop_time = amrex::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

        amrex::Print() << "Run time = " << stop_time << "s" << std::endl;
    }
    amrex::Finalize();
}
