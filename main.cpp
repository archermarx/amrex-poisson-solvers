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

    // Make BoxArray and Geometry
    amrex::BoxArray ba;
    amrex::Geometry geom;

    // Define lower and upper indices
    amrex::IntVect dom_lo(0,0,0);
    amrex::IntVect dom_hi(n_cell-1, n_cell-1, n_cell-1);

    // Make a single box that is the entire domain
    amrex::Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up box array "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box [0,1] in each direction
    amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    
    // Periodic in all directions
    amrex::Array<int,3> is_periodic{1,1,1};

    // This defines a Geometry object
    geom.define(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    // Extract dx from the Geometry object
    amrex::GpuArray<amrex::Real, 3> dx = geom.CellSizeArray();

    // Number of ghost cells for each array
    int Nghost = 1;

    // Number of components for each array
    int Ncomp = 1;

    // How boxes are distributed among MPI processes
    amrex::DistributionMapping dm(ba);

    // We allocate two multifabs - one stroes the old state, the other the new
    amrex::MultiFab phi_old(ba, dm, Ncomp, Nghost);
    amrex::MultiFab phi_new(ba, dm, Ncomp, Nghost);
    
    // time - starting time in the simulation
    amrex::Real time = 0.0;

    //************************************
    // Initialize data loop
    //************************************
    
    // Loop over boxes
    for (amrex::MFIter mfi(phi_old); mfi.isValid(); ++mfi) {
      const amrex::Box &bx = mfi.validbox();
      const amrex::Array4<amrex::Real>& phiOld = phi_old.array(mfi);
      
      // Set phi = 1 + e^(-(r - 0.5)^2)
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real x = (i + 0.5) * dx[0];
        amrex::Real y = (j + 0.5) * dx[1];
        amrex::Real z = (k + 0.5) * dx[2];
        amrex::Real rsquared = ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))/0.01;
        phiOld(i,j,k) = 1.0 + std::exp(-rsquared);
      });
    }

    //************************************
    // Write initial plotfile
    //************************************
    if (plot_int > 0) {
      int step = 0;
      const std::string &pltfile = amrex::Concatenate("output/plt", step, 5);
      WriteSingleLevelPlotfile(pltfile, phi_old, {"phi"}, geom, time, 0);
    }

    //************************************
    // Main time evolution loop
    //************************************
    for (int step = 1; step <= nsteps; ++step) {
      // fill periodic ghost cells
      phi_old.FillBoundary(geom.periodicity());
      
      // new_phi = old_phi + dt * Laplacian(old_phi)
      for (amrex::MFIter mfi(phi_old); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::Array4<amrex::Real>& phiOld = phi_old.array(mfi);
        const amrex::Array4<amrex::Real>& phiNew = phi_new.array(mfi);

        // Advance the data by dt
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
          
            phiNew(i,j,k) = phiOld(i,j,k) + dt * (
                (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
              + (phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
              + (phiOld(i,j,k+1) - 2.*phiOld(i,j,k) + phiOld(i,j,k-1)) / (dx[2]*dx[2])
            );
        });
      }
      // Update time
      time += dt;

      // Copy new solution into old solution
      amrex::MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);

      // Tell the I/O processor to write out what step we're doing
      amrex::Print() << "Advanced step " << step << "\n";

      // Write plotfile of the current data
      if (plot_int > 0 && step%plot_int == 0) {
        const std::string& pltfile = amrex::Concatenate("output/plt", step, 5);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, step);
      }
    }
  }
  amrex::Finalize();
  return 0;
}
