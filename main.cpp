#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

// Inputs
#include <AMReX_ParmParse.H>

// Diagnostics
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char *argv[]) {
  Initialize(argc, argv);
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
    Real dt;

    //************************************
    // Read parameter values from input data
    //************************************
    {
      // ParmParse is a way of reading inputs from the inputs file
      // pp.get means we require the inputs file to have it
      // pp.query measn we optionally need the inputs file to have it - but we must supply a default
      ParmParse pp;

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

    Print() << "Read inputs:\n";
    Print() << "    n_cell = " << n_cell << "\n";
    Print() << "    max_grid_size = " << max_grid_size << "\n";
    Print() << "    nsteps = " << nsteps << "\n";
    Print() << "    plot_int = " << plot_int << "\n";
    Print() << "    dt = " << dt << "\n";

    //************************************
    // Simulation setup 
    //************************************

    // Make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;

    // Define lower and upper indices
    IntVect dom_lo(AMREX_D_DECL(0,0,0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up box array "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box [0,1] in each direction
    RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    
    // Periodic in all directions
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // Extract dx from the Geometry object
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Number of ghost cells for each array
    int Nghost = 1;

    // Number of components for each array
    int Ncomp = 1;

    // How boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // We allocate two multifabs - one stroes the old state, the other the new
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);
    
    // time - starting time in the simulation
    Real time = 0.0;

    //************************************
    // Initialize data loop
    //************************************
    
    // Loop over boxes
    for (MFIter mfi(phi_old); mfi.isValid(); ++mfi) {
      const Box &bx = mfi.validbox();
      const Array4<Real>& phiOld = phi_old.array(mfi);
      
      // Set phi = 1 + e^(-(r - 0.5)^2)
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Included in all cases
        Real x = (i + 0.5) * dx[0];
        Real y = (j + 0.5) * dx[1];

        // Dimension-dependent code
#if (AMREX_SPACEDIM == 2)
        Real rsquared = ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))/0.01;
#elif (AMREX_SPACEDIM == 3)
        Real z = (k + 0.5) * dx[2];
        Real rsquared = ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))/0.01;
#endif

        // Included in all cases
        phiOld(i,j,k) = 1.0 + std::exp(-rsquared);
      });
    }

    //************************************
    // Write initial plotfile
    //************************************
    if (plot_int > 0) {
      int step = 0;
      const std::string &pltfile = Concatenate("output/plt", step, 5);
      WriteSingleLevelPlotfile(pltfile, phi_old, {"phi"}, geom, time, 0);
    }

    //************************************
    // Main time evolution loop
    //************************************
    for (int step = 1; step <= nsteps; ++step) {
      // fill periodic ghost cells
      phi_old.FillBoundary(geom.periodicity());
      
      // new_phi = old_phi + dt * Laplacian(old_phi)
      for (MFIter mfi(phi_old); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Array4<Real>& phiOld = phi_old.array(mfi);
        const Array4<Real>& phiNew = phi_new.array(mfi);

        // Advance the data by dt
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
          
            phiNew(i,j,k) = phiOld(i,j,k) + dt * (
                (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
              + (phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
              + (phiOld(i,j,k+1) - 2.*phiOld(i,j,k) + phiOld(i,j,k-1)) / (dx[2]*dx[2])
#endif
            );
        });
      }
      // Update time
      time += dt;

      // Copy new solution into old solution
      MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);

      // Tell the I/O processor to write out what step we're doing
      Print() << "Advanced step " << step << "\n";

      // Write plotfile of the current data
      if (plot_int > 0 && step%plot_int == 0) {
        const std::string& pltfile = Concatenate("output/plt", step, 5);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, step);
      }
    }
  }
  Finalize();
  return 0;
}
