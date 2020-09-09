#include "atlas/library/Library.h"
#include "atlas/projection.h"
#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"
#include "atlas/functionspace.h"

int main (int argc, char * argv[])
{
  atlas::Library::instance ().initialise (argc, argv);

  auto & comm = atlas::mpi::comm ();
  const int myproc = comm.rank ();


  using namespace atlas::util;
  using namespace atlas;
  using namespace atlas::grid;

  const int Nx = 60, Ny = 60;
  const double xmin = -20, xmax = +20, ymin = 20, ymax = +60;

  std::vector<Spacing> spacings (Ny);

  for (int i = 0; i < Ny; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", Nx) 
                         | Config ("start", xmin) | Config ("end", xmax));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "linear") | Config ("N", Ny) | Config ("start", ymin) | Config ("end", ymax));
  Projection proj (Config ("type", "lonlat"));

  atlas::StructuredGrid grid (xspace, yspace, proj, Domain ());
  atlas::grid::Distribution dist (grid, atlas::util::Config ("type", "checkerboard"));
  atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 5) | atlas::util::Config ("periodic_points", false));

  printf (" myproc = %8d\n", myproc);

  printf (" fs.j_begin () = %8d,  fs.j_end () = %8d\n", fs.j_begin (), fs.j_end ());
  for (int j = fs.j_begin (); j < fs.j_end (); j++)
    printf (" %8d > %8d .. %8d [%8d]\n", j, fs.i_begin (j), fs.i_end (j), grid.nx (j));
  printf ("-----\n");


  printf (" fs.j_begin_halo () = %8d,  fs.j_end_halo () = %8d\n", fs.j_begin_halo (), fs.j_end_halo ());
  for (int j = fs.j_begin_halo (); j < fs.j_end_halo (); j++)
    printf (" %8d > %8d .. %8d [%8d]\n", j, fs.i_begin_halo (j), fs.i_end_halo (j), grid.nx (j));


  auto f = atlas::Field ("field",
                         atlas::array::DataType::kind<double> (), 
                         atlas::array::make_shape (fs.size ()));

  
  auto v = atlas::array::make_view<double,1> (f);

  for (int i = 0; i < f.size (); i++)
    v (i) = static_cast<double> (myproc);

  fs.haloExchange (f);

  printf ("-----\n");
  for (int j = fs.j_begin_halo (); j < fs.j_end_halo (); j++)
    {
      printf (" %8d > ", j);
      for (int i = fs.i_begin_halo (j); i < fs.i_end_halo (j); i++)
        printf (" %8d", static_cast<int> (v (fs.index (i, j))));
      printf ("\n");
    }
   

  atlas::Library::instance ().finalise ();

  return 0;
}

