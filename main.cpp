#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

// Inputs
#include <AMReX_ParmParse.H>

// Diagnostics
#include <AMReX_PlotFileUtil.H>


int main (int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  {
    //************************************
    // Declare simulation parameters
    //************************************

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;
    
    // timestep
    amrex::Real dt;

    //************************************
    // Read parameter values from input data
    //************************************
    {
      // ParmParse is a way of reading inputs from the inputs file
      // pp.get means we require the inputs file to have it
      // pp.query measn we optionally need the inputs file to have it - but we must supply a default
      amrex::ParmParse pp;

      // We need to get n_cell from the inputs file - this is the number of cells on each side of a square/cubic domain
      pp.get("n_cell", n_cell);
      
      // The domain is broken up into boxes of size max_grid_size
      pp.get("max_grid_size", max_grid_size);

      // Default nsteps to 10, allow us to set it to something else in the inputs file
      nsteps = 10;
      pp.query("nsteps", nsteps);
      
      // Default plot_int to -1, allow us to set it to something else in the inputs file
      // If plot_int < 0 then no plot files will be written
      plot_int = -1;
      pp.query("plot_int", plot_int);

      // time step
      pp.get("dt", dt);
    }

    amrex::Print() << "Read inputs:\n";
    amrex::Print() << "    n_cell = " << n_cell << "\n";
    amrex::Print() << "    max_grid_size = " << max_grid_size << "\n";
    amrex::Print() << "    nsteps = " << nsteps << "\n";
    amrex::Print() << "    plot_int = " << plot_int << "\n";
    amrex::Print() << "    dt = " << dt << "\n";

    //************************************
    // Simulation setup 
    //************************************
  }
  amrex::Finalize();
  return 0;
}
