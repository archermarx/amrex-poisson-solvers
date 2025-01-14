#include "PoissonProblem.H"
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void
PoissonProblem::write_plotfile () const
{
    const int ncomp = 4;
    Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error"};

    const int nlevels = max_level+1;

    Vector<MultiFab> plotmf(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
        amrex::average_node_to_cellcenter(plotmf[ilev], 0, solution      [ilev], 0, 1);
        amrex::average_node_to_cellcenter(plotmf[ilev], 1, rhs           [ilev], 0, 1);
        amrex::average_node_to_cellcenter(plotmf[ilev], 2, exact_solution[ilev], 0, 1);

        MultiFab error(rhs[ilev].boxArray(), rhs[ilev].DistributionMap(), 1, 0);
        MultiFab::Copy(error, solution[ilev], 0, 0, 1, 0);
        MultiFab::Subtract(error, exact_solution[ilev], 0, 0, 1, 0); // error = soln - exact
        amrex::average_node_to_cellcenter(plotmf[ilev], 3, error, 0, 1);
    }

    WriteMultiLevelPlotfile("plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                            varname, geom, 0.0, Vector<int>(nlevels, 0),
                            Vector<IntVect>(nlevels, IntVect{ref_ratio}));
}
