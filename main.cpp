#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

// Diagnostics
#include <AMReX_PlotFileUtil.H>

/*
 * MultiFabs
 * C++ class for operating on and storing on multidimensional data
 * "fancy 2D or 3D vector/matrix"
 * AMReX can automatically parallize operations on data in the multifab -> easy syntax
 * from https://www.youtube.com/watch?v=498VdW2cNB8 and https://amrex-codes.github.io/amrex/tutorials_html/MultiFab.html#multifab-tutorial
 */

int main (int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  {
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

    // define a MultiFab
    // Need four elements - 
    //  1. ba - a BoxArray - list of boxes that cover our domain
    //  2. dm - a DistributionMapping - maps each box to an MPI process
    //  3. ncomp - number of components - number of data components to store in MultiFab (i.e if storing temperature and pressure, would be 2)
    //  4. ngrow - number of layers of ghost cells around the boundary (depends on application, use zero here)

    int ncomp = 1;
    int ngrow = 0;
    int ncell = 32;
    int max_grid_size = 16;
    
    // Define box array - need to define some aspects of the domain
    amrex::IntVect dom_lo(0,0,0);
    amrex::IntVect dom_hi(ncell-1, ncell-1, ncell-1);

    amrex::Box domain(dom_lo, dom_hi);
    amrex::BoxArray ba(domain);
    
    // Define maximum number of cells in each direction per MPI process
    ba.maxSize(max_grid_size);

    // Use default distribution mapping
    amrex::DistributionMapping dm(ba);

    // Create the multifab
    amrex::MultiFab mf(ba, dm, ncomp, ngrow);

    // Add data to the multifab
    amrex::RealBox real_box({0., 0., 0.}, {1., 1., 1.});
    amrex::Geometry geom(domain, &real_box);
    amrex::GpuArray<amrex::Real, 3> dx = geom.CellSizeArray();
    
    // Iterate over multifab
    // Dispatches to each MPI process, only does work on boxes that belong to a given MPI process
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      // find out local box
      const amrex::Box& bx = mfi.validbox();
      const amrex::Array4<amrex::Real>& mf_array = mf.array(mfi);

      // parallel for loop over cells/nodes in multifab
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
          
          // Get cell center coordinates
          amrex::Real x = (i + 0.5) * dx[0];
          amrex::Real y = (j + 0.5) * dx[1];
          amrex::Real z = (k + 0.5) * dx[2];
          
          // Fill with data
          amrex::Real dist_x = x - 0.5;
          amrex::Real dist_y = y - 0.5;
          amrex::Real dist_z = z - 0.5;
          amrex::Real rsquared = (dist_x*dist_x + dist_y*dist_y + dist_z*dist_z) / 0.01;            
          mf_array(i, j, k) = 1.0 + std::exp(-rsquared);
      });
    }

    // plot the MultiFab
    WriteSingleLevelPlotfile("plt001", mf, {"comp0"}, geom, 0., 0);

  }
  amrex::Finalize();
}
